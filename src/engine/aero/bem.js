/**
 * Steady blade-element-momentum (BEM) rotor aerodynamics and actuator-disc helpers.
 *
 * Formulation: one-dimensional residual in the inflow angle phi with bracketed root solve (guaranteed convergence),
 * Ning, S.A. (2014), "A simple solution method for the blade element momentum equations with guaranteed
 * convergence", Wind Energy 17:1327-1345. Prandtl tip and hub loss; high-thrust (a > 0.4) closure of Buhl, M.L.
 * (2005), "A new empirical relationship between thrust coefficient and induction factor for the turbulent windmill
 * state", NREL/TP-500-36834:  CT = 8/9 + (4F - 40/9) a + (50/9 - 4F) a^2.
 * The same code serves tidal rotors: pass water density and viscosity.
 *
 * Geometry conventions: blade stations are measured along the blade; precone and tilt are magnitudes of an upwind
 * rotor; azimuth psi = 0 with the blade pointing up, positive in the direction of rotation. Thrust is along the
 * shaft, torque about the shaft. Coefficients are referred to the swept area pi (R cos(precone))^2 and the
 * undisturbed hub-height speed U; tip-speed ratio is Omega R cos(precone) / U.
 */

const D2R = Math.PI / 180;
const TWO_PI = 2 * Math.PI;

function need(cond, name, msg) {
  if (!cond) throw new RangeError(`${name}: ${msg}`);
}

/** Betz limit of the power coefficient of an ideal actuator disc, 16/27 [-] (Betz 1920). */
export const BETZ_LIMIT = 16 / 27;

/**
 * One-dimensional momentum (actuator-disc) theory: CT = 4 a (1 - a), CP = 4 a (1 - a)^2, far-wake speed
 * U (1 - 2a). Valid for 0 <= a <= 0.5 (the wake reverses beyond); CP peaks at 16/27 for a = 1/3.
 * Reference: Burton et al., Wind Energy Handbook, section 3.2.
 * @param {number} a axial induction factor [-]
 * @returns {{ct:number, cp:number, wakeSpeedRatio:number}} thrust and power coefficients [-], far-wake/free speed [-]
 */
export function actuatorDisc(a) {
  need(a >= 0 && a <= 0.5, 'a', 'must be in [0, 0.5]');
  return { ct: 4 * a * (1 - a), cp: 4 * a * (1 - a) * (1 - a), wakeSpeedRatio: 1 - 2 * a };
}

/**
 * Thrust coefficient versus axial induction including the Buhl (2005) high-induction closure:
 * CT = 4 a F (1 - a) for a <= 0.4, CT = 8/9 + (4F - 40/9) a + (50/9 - 4F) a^2 for a > 0.4 (continuous in value and
 * slope at a = 0.4 for any F; CT = 2 at a = 1).
 * @param {number} a axial induction factor [-], 0..1
 * @param {number} [F=1] Prandtl loss factor [-], (0, 1]
 * @returns {number} local thrust coefficient [-]
 */
export function thrustCoefficientBuhl(a, F = 1) {
  need(a >= 0 && a <= 1, 'a', 'must be in [0, 1]');
  need(F > 0 && F <= 1, 'F', 'must be in (0, 1]');
  return a <= 0.4 ? 4 * a * F * (1 - a) : 8 / 9 + (4 * F - 40 / 9) * a + (50 / 9 - 4 * F) * a * a;
}

/**
 * Inverse of thrustCoefficientBuhl: axial induction for a given local thrust coefficient.
 * @param {number} ct thrust coefficient [-], 0..2
 * @param {number} [F=1] loss factor [-]
 * @returns {number} a [-]
 */
export function inductionFromCt(ct, F = 1) {
  need(ct >= 0 && ct <= 2, 'ct', 'must be in [0, 2]');
  if (ct <= 0.96 * F) return 0.5 * (1 - Math.sqrt(1 - ct / F));
  const A = 50 / 9 - 4 * F;
  const B = 4 * F - 40 / 9;
  const C = 8 / 9 - ct;
  if (Math.abs(A) < 1e-12) return -C / B;
  return (-B + Math.sqrt(Math.max(0, B * B - 4 * A * C))) / (2 * A);
}

/**
 * Maximum power coefficient of an ideal rotor with wake rotation (no drag, infinite blades) at tip-speed ratio
 * lambda (Glauert 1935):
 *   CP,max = (24 / lambda^2) * integral_{0.25}^{a2} [ (1 - a)(1 - 2a)(1 - 4a) / (1 - 3a) ]^2 da,
 * where a2 solves lambda^2 = (1 - a2)(1 - 4 a2)^2 / (1 - 3 a2). Tends to 16/27 for large lambda.
 * Reference: Manwell, McGowan and Rogers, Wind Energy Explained, 2nd ed., Eq. 3.90 and Table 3.2
 * (lambda = 1: 0.416, 2: 0.511, 5: 0.570, 10: 0.585).
 * @param {number} tsr tip-speed ratio [-] (> 0)
 * @returns {number} CP,max [-]
 */
export function glauertIdealCp(tsr) {
  need(tsr > 0, 'tsr', 'must be > 0');
  const g = (a) => (1 - a) * (1 - 4 * a) ** 2 / (1 - 3 * a) - tsr * tsr;
  let lo = 0.25;
  let hi = 1 / 3 - 1e-15;
  for (let i = 0; i < 200; i++) {
    const mid = 0.5 * (lo + hi);
    if (g(mid) < 0) lo = mid; else hi = mid;
  }
  const a2 = 0.5 * (lo + hi);
  // substitute x = 1 - 3a to keep the integrand polynomial/rational and integrate with composite Gauss-Legendre in ln x
  const f = (a) => ((1 - a) * (1 - 2 * a) * (1 - 4 * a) / (1 - 3 * a)) ** 2;
  const x1 = 1 - 3 * a2;
  const u0 = Math.log(x1);
  const u1 = Math.log(0.25);
  const n = 400;
  const gp = [-0.7745966692414834, 0, 0.7745966692414834];
  const gw = [5 / 9, 8 / 9, 5 / 9];
  let s = 0;
  for (let i = 0; i < n; i++) {
    const ua = u0 + (u1 - u0) * i / n;
    const ub = u0 + (u1 - u0) * (i + 1) / n;
    for (let k = 0; k < 3; k++) {
      const u = 0.5 * (ua + ub) + 0.5 * (ub - ua) * gp[k];
      const x = Math.exp(u);
      s += gw[k] * 0.5 * (ub - ua) * f((1 - x) / 3) * x / 3;
    }
  }
  return 24 / (tsr * tsr) * s;
}

/**
 * Prandtl loss factor F = (2/pi) acos(exp(-f)), f = (B/2) (R - r) / (r |sin phi|) at the tip and
 * f = (B/2) (r - Rhub) / (Rhub |sin phi|) at the hub; the product is returned.
 * Reference: Glauert (1935); Burton et al., Wind Energy Handbook, section 3.9. F -> 1 far from the blade ends and
 * for many blades, F -> 0 at r = R and r = Rhub.
 * @param {object} o
 * @param {number} o.blades number of blades [-]
 * @param {number} o.r local radius [m]
 * @param {number} o.tipRadius [m]
 * @param {number} [o.hubRadius=0] [m]; 0 disables the hub loss
 * @param {number} o.phi inflow angle [rad]
 * @returns {{F:number, Ftip:number, Fhub:number}} loss factors [-]
 */
export function prandtlLoss(o) {
  need(o && o.r > 0 && o.r <= o.tipRadius, 'r', 'must be in (0, tipRadius]');
  const s = Math.max(Math.abs(Math.sin(o.phi)), 1e-9);
  const Ftip = (2 / Math.PI) * Math.acos(Math.min(1, Math.exp(-0.5 * o.blades * (o.tipRadius - o.r) / (o.r * s))));
  let Fhub = 1;
  if (o.hubRadius > 0) {
    need(o.r >= o.hubRadius, 'r', 'must be >= hubRadius');
    Fhub = (2 / Math.PI) * Math.acos(Math.min(1, Math.exp(-0.5 * o.blades * (o.r - o.hubRadius) / (o.hubRadius * s))));
  }
  return { F: Ftip * Fhub, Ftip, Fhub };
}

/**
 * Axial speed-up factor u/U upstream of a tower from the potential flow about a circular cylinder (doublet in a
 * uniform stream): u/U = 1 + a^2 (y^2 - x^2) / (x^2 + y^2)^2, a = tower radius. On the upstream stagnation line
 * (y = 0) this is 1 - (a/x)^2. Only the streamwise component is returned. Valid outside the cylinder for upwind
 * rotors; it does not represent a downwind (wake) shadow.
 * Reference: Burton et al., Wind Energy Handbook, section 5.7.2.
 * @param {number} x streamwise distance from the tower axis [m] (either sign)
 * @param {number} y lateral distance from the tower axis [m]
 * @param {number} towerRadius local tower radius [m]
 * @returns {number} u/U [-]
 */
export function towerPotentialFactor(x, y, towerRadius) {
  const d2 = x * x + y * y;
  need(towerRadius >= 0, 'towerRadius', 'must be >= 0');
  if (d2 <= towerRadius * towerRadius) return 0;
  return 1 + towerRadius * towerRadius * (y * y - x * x) / (d2 * d2);
}

// ---- element solve (module-level scratch: no allocation in the root search) -------------------------------------
const E = {
  Vx: 0, Vy: 0, r: 1, chord: 1, theta: 0, polar: null, B: 3, Rtip: 1, Rhub: 0,
  tipLoss: true, hubLoss: true, useCd: true, wakeRot: true, buhl: true,
  a: 0, ap: 0, F: 1, cl: 0, cd: 0, alpha: 0,
};

function resid(phi) {
  const sphi = Math.sin(phi);
  const cphi = Math.cos(phi);
  const alpha = phi - E.theta;
  const cl = E.polar.clAt(alpha);
  const cd = E.polar.cdAt(alpha);
  E.alpha = alpha; E.cl = cl; E.cd = cd;
  let cn; let ct;
  if (E.useCd) { cn = cl * cphi + cd * sphi; ct = cl * sphi - cd * cphi; } else { cn = cl * cphi; ct = cl * sphi; }
  let F = 1;
  const as = Math.abs(sphi);
  if (E.tipLoss) {
    const f = 0.5 * E.B * (E.Rtip - E.r) / (E.r * as);
    F *= (2 / Math.PI) * Math.acos(Math.min(1, Math.exp(-f)));
  }
  if (E.hubLoss && E.Rhub > 0) {
    const f = 0.5 * E.B * (E.r - E.Rhub) / (E.Rhub * as);
    F *= (2 / Math.PI) * Math.acos(Math.min(1, Math.exp(-f)));
  }
  if (F < 1e-6) F = 1e-6;
  E.F = F;
  const sig = E.B * E.chord / (TWO_PI * E.r);
  const k = sig * cn / (4 * F * sphi * sphi);
  let kp = sig * ct / (4 * F * sphi * cphi);
  let a;
  if (phi > 0) {
    if (k <= 2 / 3 || !E.buhl) {
      a = k / (1 + k);
    } else {
      const g1 = 2 * F * k - (10 / 9 - F);
      const g2 = 2 * F * k - (4 / 3 - F) * F;
      const g3 = 2 * F * k - (25 / 9 - 2 * F);
      a = Math.abs(g3) < 1e-6 ? 1 - 1 / (2 * Math.sqrt(g2)) : (g1 - Math.sqrt(g2)) / g3;
    }
  } else {
    a = k > 1 ? k / (k - 1) : 0;
  }
  let ap = kp / (1 - kp);
  if (!E.wakeRot) { ap = 0; kp = 0; }
  E.a = a; E.ap = ap;
  const lr = E.Vy / E.Vx;
  return phi > 0 ? sphi / (1 - a) - cphi / lr * (1 - kp) : sphi * (1 - k) - cphi / lr * (1 - kp);
}

function brentPhi(lo, hi, flo, fhi) {
  let a = lo; let b = hi; let fa = flo; let fb = fhi;
  let c = a; let fc = fa; let d = b - a; let e = d;
  for (let it = 0; it < 100; it++) {
    if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) { c = a; fc = fa; d = b - a; e = d; }
    if (Math.abs(fc) < Math.abs(fb)) { a = b; b = c; c = a; fa = fb; fb = fc; fc = fa; }
    const tol = 2 * 2.2e-16 * Math.abs(b) + 0.5e-10;
    const xm = 0.5 * (c - b);
    if (Math.abs(xm) <= tol || fb === 0) return b;
    if (Math.abs(e) >= tol && Math.abs(fa) > Math.abs(fb)) {
      const s = fb / fa;
      let p; let q;
      if (a === c) { p = 2 * xm * s; q = 1 - s; } else {
        const qq = fa / fc; const r = fb / fc;
        p = s * (2 * xm * qq * (qq - r) - (b - a) * (r - 1));
        q = (qq - 1) * (r - 1) * (s - 1);
      }
      if (p > 0) q = -q;
      p = Math.abs(p);
      if (2 * p < Math.min(3 * xm * q - Math.abs(tol * q), Math.abs(e * q))) { e = d; d = p / q; } else { d = xm; e = d; }
    } else { d = xm; e = d; }
    a = b; fa = fb;
    b += Math.abs(d) > tol ? d : (xm > 0 ? tol : -tol);
    fb = resid(b);
  }
  return b;
}

// Solves the element with E.* set; returns phi and leaves a, ap, F, cl, cd, alpha in E. `converged` false -> no
// induction (parked / reversed inflow).
let E_converged = true;
function solvePhi() {
  E_converged = true;
  if (!(E.Vx > 1e-9) || !(E.Vy > 1e-6 * E.Vx)) {
    E_converged = false;
    const phi = Math.atan2(E.Vx, E.Vy);
    E.alpha = phi - E.theta; E.cl = E.polar.clAt(E.alpha); E.cd = E.polar.cdAt(E.alpha);
    E.a = 0; E.ap = 0; E.F = 1;
    return phi;
  }
  const eps = 1e-6;
  let lo = eps; let hi = Math.PI / 2;
  let flo = resid(lo); let fhi = resid(hi);
  if (flo * fhi > 0) {
    lo = -Math.PI / 4; hi = -eps;
    flo = resid(lo); fhi = resid(hi);
    if (flo * fhi > 0) {
      lo = Math.PI / 2; hi = Math.PI - eps;
      flo = resid(lo); fhi = resid(hi);
      if (flo * fhi > 0) {
        E_converged = false;
        const phi = Math.atan2(E.Vx, E.Vy);
        E.alpha = phi - E.theta; E.cl = E.polar.clAt(E.alpha); E.cd = E.polar.cdAt(E.alpha);
        E.a = 0; E.ap = 0; E.F = 1;
        return phi;
      }
    }
  }
  const phi = brentPhi(lo, hi, flo, fhi);
  resid(phi);
  return phi;
}

/**
 * Solve one blade element (annulus) for the induction factors.
 * Residual of Ning (2014): R(phi) = sin(phi)/(1 - a) - (Vx/Vy) cos(phi) (1 - k'), with
 * k = sigma' cn / (4 F sin^2 phi), k' = sigma' ct / (4 F sin phi cos phi), a = k/(1 + k) (momentum, k <= 2/3) or the
 * Buhl closure beyond, a' = k'/(1 - k'); brackets (0, pi/2], [-pi/4, 0) and (pi/2, pi) are tried in that order.
 * When the element has no rotation-induced inflow (Vy <= 0) or no bracket exists the element is evaluated without
 * induction and `converged` is false.
 * @param {object} o
 * @param {number} o.Vx inflow normal to the rotor plane at the element, without induction [m/s]
 * @param {number} o.Vy inflow tangential to the rotor plane (Omega r plus in-plane wind) [m/s]
 * @param {number} o.r local radius [m]
 * @param {number} o.chord chord [m]
 * @param {number} o.theta twist + pitch, angle from rotor plane to chord line [rad]
 * @param {object} o.polar polar with clAt/cdAt (polar.js)
 * @param {number} o.blades number of blades [-]
 * @param {number} o.tipRadius [m]
 * @param {number} [o.hubRadius=0] [m]
 * @param {number} [o.rho=1.225] fluid density [kg/m3]
 * @param {boolean} [o.tipLoss=true]
 * @param {boolean} [o.hubLoss=true]
 * @param {boolean} [o.dragInInduction=true] include drag in the momentum balance
 * @param {boolean} [o.wakeRotation=true] solve tangential induction
 * @param {boolean} [o.highThrust=true] Buhl closure for a > 0.4 (false: plain momentum theory at all loadings)
 * @param {object} [out] optional object to receive the result fields (avoids allocation in time marching)
 * @returns {{phi:number,a:number,ap:number,alpha:number,cl:number,cd:number,cn:number,ct:number,F:number,W:number,
 *   Np:number,Tp:number,converged:boolean}} inflow angle [rad], induction factors [-], angle of attack [rad],
 *   coefficients [-], loss factor [-], relative speed [m/s], normal and tangential force per unit length [N/m]
 */
export function bemElement(o, out = {}) {
  E.Vx = o.Vx; E.Vy = o.Vy; E.r = o.r; E.chord = o.chord; E.theta = o.theta; E.polar = o.polar; E.B = o.blades;
  E.Rtip = o.tipRadius; E.Rhub = o.hubRadius ?? 0;
  E.tipLoss = o.tipLoss ?? true; E.hubLoss = o.hubLoss ?? true; E.useCd = o.dragInInduction ?? true;
  E.wakeRot = o.wakeRotation ?? true; E.buhl = o.highThrust ?? true;
  const phi = solvePhi();
  const s = Math.sin(phi); const c = Math.cos(phi);
  const W = relSpeed(phi);
  const q = 0.5 * (o.rho ?? 1.225) * W * W * o.chord;
  out.phi = phi; out.a = E.a; out.ap = E.ap; out.alpha = E.alpha; out.cl = E.cl; out.cd = E.cd;
  out.cn = E.cl * c + E.cd * s; out.ct = E.cl * s - E.cd * c; out.F = E.F; out.W = W;
  out.Np = out.cn * q; out.Tp = out.ct * q; out.converged = E_converged;
  return out;
}

function relSpeed(phi) {
  const vx = E.Vx * (1 - E.a);
  const vy = E.Vy * (1 + E.ap);
  if (Math.abs(E.a) > 10) return Math.abs(vy / Math.cos(phi));
  if (Math.abs(E.ap) > 10) return Math.abs(vx / Math.sin(phi));
  return Math.sqrt(vx * vx + vy * vy);
}

function towerRadiusAt(ts, zBelowTop) {
  if (zBelowTop < 0) return 0;
  const s = Math.min(1, zBelowTop / ts.height);
  return 0.5 * (ts.topDiameter + (ts.baseDiameter - ts.topDiameter) * s);
}

/**
 * Steady BEM solution of a rotor.
 *
 * Inflow at a blade element at azimuth psi (Ning 2014 / CCBlade kinematics, straight blades):
 *   Vx = V [ (cos(yaw) sin(tilt) cos(psi) + sin(yaw) sin(psi)) sin(cone) + cos(yaw) cos(tilt) cos(cone) ]
 *   Vy = V [ cos(yaw) sin(tilt) sin(psi) - sin(yaw) cos(psi) ] + Omega r cos(cone)
 * with V the sheared wind at the element height, V = U (1 + z/hubHeight)^shearExp, times the tower potential-flow
 * factor when `towerShadow` is given. Loads are averaged over `nSectors` azimuth positions.
 * Skewed wake (yaw and/or tilt): after the axisymmetric solve the axial induction is redistributed as
 *   a_skew = a [1 + (15 pi / 32) (r/R) tan(chi/2) cos(psi_s)],  chi = (0.6 a + 1) gamma,
 * gamma = angle between wind and shaft, psi_s measured from the most downwind blade position (Pitt-Peters / Glauert
 * form with the wake-skew approximation of Burton et al., as used in AeroDyn: Moriarty and Hansen (2005),
 * NREL/TP-500-36881); reliable for gamma below about 30 deg.
 * Integration along the blade is trapezoidal through the stations with zero loading at hub and tip.
 * Limits: steady, rigid, straight blades, no radial flow; deep-stall and propeller-brake states are evaluated but
 * are outside the validity of momentum theory.
 *
 * Parked and idling rotor. Momentum theory presupposes a rotor that extracts energy from a stream tube; a parked
 * or slowly idling rotor (blades feathered, tip-speed ratio near zero) does not form one, and the induction
 * iteration then returns loads that jump between solution branches as the inflow angle changes with yaw. Below
 * the tip-speed ratio `parkedTsr` (default 0.2, below the lowest tip-speed ratio of the default cpSurface grid) the
 * induction is therefore switched off: each element sees the
 * geometric inflow (wind components normal and tangential to the rotor plane at its azimuth, plus Omega r), the
 * lift and drag of its polar at that angle of attack give the forces, and the loads are averaged over at least 36
 * azimuth positions — the mean over rotor positions of a rotor stopped at an arbitrary azimuth; `azimuthal`
 * gives the single-blade loads per position. Tip and hub losses are not applied (F = 1) and the result carries
 * `parked: true`. This is quasi-steady strip theory on stationary blades: no dynamic stall, no spanwise flow, and
 * the polars must cover the full ±180 deg range. `parkedTsr: 0` restores the momentum solution everywhere.
 * @param {object} o
 * @param {object} o.rotor rotor (rotors.js createRotor)
 * @param {number} o.U undisturbed hub-height speed [m/s] (>= 0)
 * @param {number} o.omega rotor speed [rad/s] (>= 0)
 * @param {number} [o.pitchDeg=0] collective pitch, positive to feather [deg]
 * @param {number} [o.rho=1.225] density [kg/m3] (sea water 1025)
 * @param {number} [o.mu=1.81e-5] dynamic viscosity [Pa s] (sea water about 1.4e-3), used for the reported Reynolds number
 * @param {number} [o.yawDeg=0] yaw misalignment [deg]
 * @param {number} [o.tiltDeg] shaft tilt [deg]; default rotor.tiltDeg
 * @param {number} [o.preconeDeg] precone [deg]; default rotor.preconeDeg
 * @param {number} [o.shearExp=0] power-law shear exponent [-]
 * @param {number} [o.hubHeight] [m]; default rotor.hubHeight (needed when shearExp != 0)
 * @param {object|boolean} [o.towerShadow=false] true (use rotor.tower and rotor.overhang) or
 *   { topDiameter [m], baseDiameter [m], height [m], overhang [m], hubAboveTowerTop [m] = 0 }
 * @param {number} [o.nSectors] azimuth sectors; default 1 (axisymmetric), 4 (tilt only) or 12 (shear, yaw, shadow)
 * @param {boolean} [o.tipLoss=true]
 * @param {boolean} [o.hubLoss=true]
 * @param {boolean} [o.dragInInduction=true]
 * @param {boolean} [o.wakeRotation=true]
 * @param {boolean} [o.highThrust=true]
 * @param {boolean} [o.skewCorrection=true]
 * @param {number} [o.parkedTsr=0.2] tip-speed ratio below which the parked-rotor branch (no induction) is used [-]
 * @returns {{thrust:number,torque:number,power:number,Cp:number,Ct:number,Cq:number,tsr:number,parked:boolean,
 *   flapMoment:number,edgeMoment:number,flapMomentMax:number,converged:boolean,
 *   spanwise:Array<{r:number,a:number,ap:number,phi:number,alpha:number,cl:number,cd:number,Np:number,Tp:number,Re:number,F:number,W:number}>,
 *   azimuthal:{psi:number[],bladeThrust:number[],bladeTorque:number[],flapMoment:number[],edgeMoment:number[]}}}
 *   thrust [N], torque [N m], aerodynamic power torque*omega [W], coefficients [-], tip-speed ratio [-];
 *   azimuth-mean aerodynamic blade-root flapwise (out-of-plane) and edgewise (in-plane) moments about the hub radius
 *   [N m] and the largest flapwise value over the azimuth; spanwise azimuth-averaged quantities (angles in rad,
 *   forces per length N/m); per-blade loads versus azimuth [rad] (thrust N, torque N m, moments N m)
 */
export function bemSolve(o) {
  const rotor = o.rotor;
  need(rotor && Array.isArray(rotor.stations), 'rotor', 'rotor with stations required');
  need(o.U >= 0, 'U', 'must be >= 0');
  need(o.omega >= 0, 'omega', 'must be >= 0');
  const rho = o.rho ?? 1.225;
  const mu = o.mu ?? 1.81e-5;
  need(rho > 0, 'rho', 'must be > 0');
  const pitch = (o.pitchDeg ?? 0) * D2R;
  const yaw = (o.yawDeg ?? 0) * D2R;
  const tilt = (o.tiltDeg ?? rotor.tiltDeg ?? 0) * D2R;
  const cone = (o.preconeDeg ?? rotor.preconeDeg ?? 0) * D2R;
  const shear = o.shearExp ?? 0;
  const hubHt = o.hubHeight ?? rotor.hubHeight;
  if (shear !== 0) need(hubHt > 0, 'hubHeight', 'required when shearExp != 0');
  let ts = null;
  if (o.towerShadow) {
    ts = o.towerShadow === true
      ? { topDiameter: rotor.tower?.topDiameter, baseDiameter: rotor.tower?.baseDiameter, height: rotor.tower?.height, overhang: rotor.overhang, hubAboveTowerTop: 0 }
      : { hubAboveTowerTop: 0, ...o.towerShadow };
    need(ts.topDiameter > 0 && ts.baseDiameter > 0 && ts.height > 0 && ts.overhang > 0, 'towerShadow', 'topDiameter, baseDiameter, height and overhang required');
  }
  const asym = shear !== 0 || yaw !== 0 || ts !== null;
  const parkedTsr = o.parkedTsr ?? 0.2;
  need(parkedTsr >= 0, 'parkedTsr', 'must be >= 0');
  const parked = o.U > 0 && o.omega * rotor.tipRadius < parkedTsr * o.U;
  const nSec = o.nSectors ?? (parked ? 36 : (asym ? 12 : (tilt !== 0 ? 4 : 1)));
  need(Number.isInteger(nSec) && nSec >= 1, 'nSectors', 'integer >= 1 required');
  const st = rotor.stations;
  const n = st.length;
  const Rtip = rotor.tipRadius;
  const Rhub = rotor.hubRadius;
  const B = rotor.blades;
  const sy = Math.sin(yaw); const cy = Math.cos(yaw);
  const sT = Math.sin(tilt); const cT = Math.cos(tilt);
  const sc = Math.sin(cone); const cc = Math.cos(cone);
  const Rrot = Rtip * cc;
  const area = Math.PI * Rrot * Rrot;
  const skewOn = (o.skewCorrection ?? true);
  const inPlane = Math.sqrt(cy * cy * sT * sT + sy * sy);
  const gamma = Math.atan2(inPlane, cy * cT);

  E.B = B; E.Rtip = Rtip; E.Rhub = Rhub;
  E.tipLoss = o.tipLoss ?? true; E.hubLoss = o.hubLoss ?? true; E.useCd = o.dragInInduction ?? true;
  E.wakeRot = o.wakeRotation ?? true; E.buhl = o.highThrust ?? true;

  const acc = new Float64Array(n * 11);
  const Np = new Float64Array(n);
  const Tp = new Float64Array(n);
  const az = { psi: [], bladeThrust: [], bladeTorque: [], flapMoment: [], edgeMoment: [] };
  let thrust = 0; let torque = 0; let flap = 0; let edge = 0; let flapMax = -Infinity; let allConv = true;

  for (let j = 0; j < nSec; j++) {
    const psi = TWO_PI * j / nSec;
    const sa = Math.sin(psi); const ca = Math.cos(psi);
    for (let i = 0; i < n; i++) {
      const s = st[i];
      const r = s.r;
      let V = o.U;
      const z = r * cc * ca * cT + r * sc * sT;          // element height above hub [m]
      if (shear !== 0) V *= Math.pow(Math.max(1e-3, 1 + z / hubHt), shear);
      if (ts) {
        const zBelowTop = -z - ts.hubAboveTowerTop;
        if (zBelowTop > 0) {
          const xUp = ts.overhang * cT + r * (sc * cT - cc * ca * sT);
          V *= towerPotentialFactor(xUp, r * cc * sa, towerRadiusAt(ts, zBelowTop));
        }
      }
      const radial = cy * sT * ca + sy * sa;             // wind component along the blade / V
      E.Vx = V * (radial * sc + cy * cT * cc);
      E.Vy = V * (cy * sT * sa - sy * ca) + o.omega * r * cc;
      E.r = r; E.chord = s.chord; E.theta = s.twistDeg * D2R + pitch; E.polar = s.polar;
      let phi;
      if (parked) {
        // stationary blade element in the geometric inflow: no induction, no tip or hub loss
        phi = Math.atan2(E.Vx, E.Vy);
        E.a = 0; E.ap = 0; E.F = 1;
        E.alpha = phi - E.theta; E.cl = s.polar.clAt(E.alpha); E.cd = s.polar.cdAt(E.alpha);
      } else {
        phi = solvePhi();
        if (!E_converged && o.omega > 0 && o.U > 0) allConv = false;
      }
      if (!parked && skewOn && E_converged && gamma > 1e-6 && E.a > 0 && E.a < 1) {
        const chi = (0.6 * E.a + 1) * gamma;
        const fac = 1 + (15 * Math.PI / 32) * (r / Rtip) * Math.tan(0.5 * chi) * (radial / inPlane);
        E.a = Math.min(0.95, E.a * fac);
        phi = Math.atan2(E.Vx * (1 - E.a), E.Vy * (1 + E.ap));
        E.alpha = phi - E.theta; E.cl = s.polar.clAt(E.alpha); E.cd = s.polar.cdAt(E.alpha);
      }
      const sp = Math.sin(phi); const cp = Math.cos(phi);
      const W = relSpeed(phi);
      const q = 0.5 * rho * W * W * s.chord;
      Np[i] = (E.cl * cp + E.cd * sp) * q;
      Tp[i] = (E.cl * sp - E.cd * cp) * q;
      const k = i * 11;
      acc[k] += E.a; acc[k + 1] += E.ap; acc[k + 2] += phi; acc[k + 3] += E.alpha; acc[k + 4] += E.cl; acc[k + 5] += E.cd;
      acc[k + 6] += Np[i]; acc[k + 7] += Tp[i]; acc[k + 8] += rho * W * s.chord / mu; acc[k + 9] += E.F; acc[k + 10] += W;
    }
    // trapezoidal integration hub -> stations -> tip with zero end loading
    let T1 = 0; let Q1 = 0; let Mf = 0; let Me = 0;
    let rPrev = Rhub; let nPrev = 0; let tPrev = 0;
    for (let i = 0; i <= n; i++) {
      const r = i < n ? st[i].r : Rtip;
      const nn = i < n ? Np[i] : 0;
      const tt = i < n ? Tp[i] : 0;
      const dr = r - rPrev;
      T1 += 0.5 * (nn + nPrev) * dr;
      Q1 += 0.5 * (tt * r + tPrev * rPrev) * dr;
      Mf += 0.5 * (nn * (r - Rhub) + nPrev * (rPrev - Rhub)) * dr;
      Me += 0.5 * (tt * (r - Rhub) + tPrev * (rPrev - Rhub)) * dr;
      rPrev = r; nPrev = nn; tPrev = tt;
    }
    T1 *= cc; Q1 *= cc;
    az.psi.push(psi); az.bladeThrust.push(T1); az.bladeTorque.push(Q1); az.flapMoment.push(Mf); az.edgeMoment.push(Me);
    thrust += T1; torque += Q1; flap += Mf; edge += Me;
    if (Mf > flapMax) flapMax = Mf;
  }
  thrust *= B / nSec; torque *= B / nSec; flap /= nSec; edge /= nSec;
  const spanwise = new Array(n);
  for (let i = 0; i < n; i++) {
    const k = i * 11;
    spanwise[i] = {
      r: st[i].r, a: acc[k] / nSec, ap: acc[k + 1] / nSec, phi: acc[k + 2] / nSec, alpha: acc[k + 3] / nSec,
      cl: acc[k + 4] / nSec, cd: acc[k + 5] / nSec, Np: acc[k + 6] / nSec, Tp: acc[k + 7] / nSec,
      Re: acc[k + 8] / nSec, F: acc[k + 9] / nSec, W: acc[k + 10] / nSec,
    };
  }
  const power = torque * o.omega;
  const qdyn = 0.5 * rho * o.U * o.U * area;
  const has = qdyn > 0;
  return {
    thrust, torque, power,
    Cp: has ? power / (qdyn * o.U) : 0, Ct: has ? thrust / qdyn : 0, Cq: has ? torque / (qdyn * Rrot) : 0,
    tsr: has ? o.omega * Rrot / o.U : Infinity, parked,
    flapMoment: flap, edgeMoment: edge, flapMomentMax: flapMax, converged: allConv, spanwise, azimuthal: az,
  };
}

// reference radius of the coefficients for a rotor and bemSolve overrides: tip radius x cos(precone)
function refRadius(rotor, bem) {
  const cone = ((bem && bem.preconeDeg) ?? rotor.preconeDeg ?? 0) * D2R;
  return rotor.tipRadius * Math.cos(cone);
}

function range(a, b, step) {
  const out = [];
  const n = Math.round((b - a) / step);
  for (let i = 0; i <= n; i++) out.push(a + i * step);
  return out;
}

/**
 * Cp, Ct and Cq surfaces over tip-speed ratio and collective pitch (steady BEM at a reference wind speed; with
 * Reynolds-independent polars the coefficients depend on tsr and pitch only).
 * @param {object} o
 * @param {object} o.rotor rotor
 * @param {number[]} [o.tsr] tip-speed ratios [-], ascending; default 0.25, 0.5, then 1..18 step 0.5
 * @param {number[]} [o.pitchDeg] pitch angles [deg], ascending; default -5..30 step 1, then 35..90 step 5 (so that
 *   feathering and idling are covered; lookups are clamped outside the table)
 * @param {number} [o.Uref=10] wind speed at which the BEM is run [m/s]
 * @param {number} [o.rho=1.225] density [kg/m3]
 * @param {object} [o.bem] further bemSolve options (yawDeg, tiltDeg, shearExp, nSectors, tipLoss, ...)
 * @returns {{tsr:number[],pitchDeg:number[],cp:number[][],ct:number[][],cq:number[][],R:number,area:number}}
 *   coefficient tables indexed [iPitch][iTsr] [-], reference radius [m] and area [m2] of the coefficients
 */
export function cpSurface(o) {
  const tsr = o.tsr ?? [0.25, 0.5, ...range(1, 18, 0.5)];
  const pit = o.pitchDeg ?? [...range(-5, 30, 1), ...range(35, 90, 5)];
  const Uref = o.Uref ?? 10;
  need(tsr.length >= 2 && pit.length >= 2, 'tsr/pitchDeg', 'at least two values each');
  for (let i = 1; i < tsr.length; i++) need(tsr[i] > tsr[i - 1], 'tsr', 'must be ascending');
  for (let i = 1; i < pit.length; i++) need(pit[i] > pit[i - 1], 'pitchDeg', 'must be ascending');
  need(tsr[0] > 0, 'tsr', 'must be > 0');
  const R = refRadius(o.rotor, o.bem);
  const cp = []; const ct = []; const cq = [];
  const args = { ...(o.bem || {}), rotor: o.rotor, U: Uref, rho: o.rho ?? 1.225, omega: 0, pitchDeg: 0 };
  for (let ip = 0; ip < pit.length; ip++) {
    const rp = new Array(tsr.length); const rt = new Array(tsr.length); const rq = new Array(tsr.length);
    args.pitchDeg = pit[ip];
    for (let it = 0; it < tsr.length; it++) {
      args.omega = tsr[it] * Uref / R;
      const s = bemSolve(args);
      rp[it] = s.Cp; rt[it] = s.Ct; rq[it] = s.Cq;
    }
    cp.push(rp); ct.push(rt); cq.push(rq);
  }
  return { tsr: tsr.slice(), pitchDeg: pit.slice(), cp, ct, cq, R, area: Math.PI * R * R };
}

function cell(xs, x) {
  const n = xs.length;
  if (x <= xs[0]) return [0, 0];
  if (x >= xs[n - 1]) return [n - 2, 1];
  let lo = 0; let hi = n - 1;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] <= x) lo = m; else hi = m; }
  return [lo, (x - xs[lo]) / (xs[hi] - xs[lo])];
}

/**
 * Bilinear lookup in a cpSurface result (clamped at the table edges), with partial derivatives of the cell.
 * @param {object} surface cpSurface result
 * @param {number} tsr tip-speed ratio [-]
 * @param {number} pitchDeg pitch [deg]
 * @returns {{cp:number, ct:number, cq:number, dCpDTsr:number, dCpDPitchDeg:number, dCtDTsr:number, dCtDPitchDeg:number}}
 *   coefficients [-] and their slopes per unit tsr [-] and per degree of pitch [1/deg]
 */
export function surfaceLookup(surface, tsr, pitchDeg) {
  const [i, s] = cell(surface.tsr, tsr);
  const [j, t] = cell(surface.pitchDeg, pitchDeg);
  const dx = surface.tsr[i + 1] - surface.tsr[i];
  const dy = surface.pitchDeg[j + 1] - surface.pitchDeg[j];
  const bil = (z) => {
    const z00 = z[j][i]; const z10 = z[j][i + 1]; const z01 = z[j + 1][i]; const z11 = z[j + 1][i + 1];
    return [
      (1 - s) * (1 - t) * z00 + s * (1 - t) * z10 + (1 - s) * t * z01 + s * t * z11,
      ((1 - t) * (z10 - z00) + t * (z11 - z01)) / dx,
      ((1 - s) * (z01 - z00) + s * (z11 - z10)) / dy,
    ];
  };
  const p = bil(surface.cp); const c = bil(surface.ct); const q = bil(surface.cq);
  return { cp: p[0], ct: c[0], cq: q[0], dCpDTsr: p[1], dCpDPitchDeg: p[2], dCtDTsr: c[1], dCtDPitchDeg: c[2] };
}

function golden(f, a, b, tol) {
  const g = (Math.sqrt(5) - 1) / 2;
  let c = b - g * (b - a); let d = a + g * (b - a);
  let fc = f(c); let fd = f(d);
  while (b - a > tol) {
    if (fc < fd) { b = d; d = c; fd = fc; c = b - g * (b - a); fc = f(c); } else { a = c; c = d; fc = fd; d = a + g * (b - a); fd = f(d); }
  }
  return 0.5 * (a + b);
}

/**
 * Maximum power coefficient and the optimum tip-speed ratio and pitch (coarse grid search followed by alternating
 * golden-section refinement in tsr and pitch).
 * @param {object} o
 * @param {object} o.rotor rotor
 * @param {number[]} [o.tsrRange=[3, 14]] search range [-]
 * @param {number[]} [o.pitchRangeDeg=[-5, 8]] search range [deg]
 * @param {number|null} [o.pitchDeg=null] fix the pitch [deg] and optimise the tip-speed ratio only
 * @param {number} [o.Uref=10] wind speed for the BEM runs [m/s]
 * @param {number} [o.rho=1.225] [kg/m3]
 * @param {object} [o.bem] further bemSolve options
 * @returns {{cpMax:number, tsrOpt:number, pitchOptDeg:number, ct:number}} maximum Cp [-], optimum tsr [-] and pitch
 *   [deg], thrust coefficient at the optimum [-]
 */
export function maxCp(o) {
  const tr = o.tsrRange ?? [3, 14];
  const pr = o.pitchRangeDeg ?? [-5, 8];
  const Uref = o.Uref ?? 10;
  const R = refRadius(o.rotor, o.bem);
  const args = { ...(o.bem || {}), rotor: o.rotor, U: Uref, rho: o.rho ?? 1.225, omega: 0, pitchDeg: 0 };
  let last = null;
  const f = (l, p) => { args.omega = l * Uref / R; args.pitchDeg = p; last = bemSolve(args); return -last.Cp; };
  const fixed = o.pitchDeg !== undefined && o.pitchDeg !== null;
  let best = [Infinity, 0, 0];
  const nl = 23;
  const np = fixed ? 1 : 14;
  for (let i = 0; i < nl; i++) {
    const l = tr[0] + (tr[1] - tr[0]) * i / (nl - 1);
    for (let j = 0; j < np; j++) {
      const p = fixed ? o.pitchDeg : pr[0] + (pr[1] - pr[0]) * j / (np - 1);
      const v = f(l, p);
      if (v < best[0]) best = [v, l, p];
    }
  }
  let l = best[1]; let p = best[2];
  const dl = (tr[1] - tr[0]) / (nl - 1);
  const dp = fixed ? 0 : (pr[1] - pr[0]) / (np - 1);
  for (let it = 0; it < 4; it++) {
    l = golden((x) => f(x, p), Math.max(tr[0], l - dl), Math.min(tr[1], l + dl), 1e-3);
    if (!fixed) p = golden((x) => f(l, x), Math.max(pr[0], p - dp), Math.min(pr[1], p + dp), 1e-3);
  }
  const cp = -f(l, p);
  return { cpMax: cp, tsrOpt: l, pitchOptDeg: p, ct: last.Ct };
}

/**
 * Steady power curve of a variable-speed pitch-regulated rotor from direct BEM solutions.
 * Schedule: below rated the rotor follows the optimum tip-speed ratio at fine pitch, with the speed clipped to
 * [minRotorSpeed, ratedRotorSpeed]; once the electrical power at rated speed would exceed rated power the pitch is
 * increased (to feather) to hold rated power. No hysteresis, losses lumped in `drivetrainEfficiency`.
 * (control.js `operatingSchedule` does the same from pre-computed surfaces and adds peak shaving.)
 * @param {object} o
 * @param {object} o.rotor rotor with ratedPower [W] and ratedRotorSpeed [rad/s]
 * @param {number[]} [o.speeds] hub-height wind speeds [m/s]; default cutIn..cutOut in 0.5 m/s steps
 * @param {number} [o.rho=1.225] [kg/m3]
 * @param {number} [o.ratedPower] electrical rated power [W]; default rotor.ratedPower
 * @param {number} [o.efficiency] drivetrain efficiency [-]; default rotor.drivetrainEfficiency or 1
 * @param {number} [o.tsrOpt] operating tip-speed ratio [-]; default from maxCp
 * @param {number} [o.pitchFineDeg] fine pitch [deg]; default from maxCp
 * @param {number} [o.minRotorSpeed] [rad/s]; default rotor.minRotorSpeed or 0
 * @param {number} [o.ratedRotorSpeed] [rad/s]; default rotor.ratedRotorSpeed
 * @param {number} [o.cutIn] [m/s]; default rotor.cutIn or 0
 * @param {number} [o.cutOut] [m/s]; default rotor.cutOut or Infinity
 * @param {object} [o.bem] further bemSolve options
 * @returns {{U:number[],power:number[],aeroPower:number[],thrust:number[],torque:number[],omega:number[],
 *   pitchDeg:number[],cp:number[],ct:number[],tsr:number[],region:string[],ratedWind:number,cpMax:number,
 *   tsrOpt:number,pitchFineDeg:number}} electrical power [W], aerodynamic power [W], thrust [N], aerodynamic torque
 *   [N m], rotor speed [rad/s], pitch [deg], coefficients [-], region label ('off', '1.5', '2', '2.5', '3'),
 *   rated wind speed [m/s] (NaN if rated power is not reached within the speed range)
 */
export function powerCurve(o) {
  const rotor = o.rotor;
  const rho = o.rho ?? 1.225;
  const Prated = o.ratedPower ?? rotor.ratedPower;
  const eta = o.efficiency ?? rotor.drivetrainEfficiency ?? 1;
  const wRated = o.ratedRotorSpeed ?? rotor.ratedRotorSpeed;
  const wMin = o.minRotorSpeed ?? rotor.minRotorSpeed ?? 0;
  const cutIn = o.cutIn ?? rotor.cutIn ?? 0;
  const cutOut = o.cutOut ?? rotor.cutOut ?? Infinity;
  need(Prated > 0, 'ratedPower', 'must be > 0');
  need(wRated > 0, 'ratedRotorSpeed', 'must be > 0');
  need(eta > 0 && eta <= 1, 'efficiency', 'must be in (0, 1]');
  let tsrOpt = o.tsrOpt; let fine = o.pitchFineDeg; let cpMax;
  if (tsrOpt === undefined || fine === undefined) {
    const m = maxCp({ rotor, rho, bem: o.bem, pitchDeg: fine });
    tsrOpt = tsrOpt ?? m.tsrOpt; fine = fine ?? m.pitchOptDeg; cpMax = m.cpMax;
  }
  const R = refRadius(rotor, o.bem);
  const speeds = o.speeds ?? range(Math.max(cutIn, 0.5), Number.isFinite(cutOut) ? cutOut : 25, 0.5);
  const args = { ...(o.bem || {}), rotor, U: 0, rho, omega: 0, pitchDeg: fine };
  const run = (U, w, p) => { args.U = U; args.omega = w; args.pitchDeg = p; return bemSolve(args); };
  const wAt = (U) => Math.min(wRated, Math.max(wMin, tsrOpt * U / R));
  if (cpMax === undefined) cpMax = run(10, tsrOpt * 10 / R, fine).Cp;
  const out = { U: [], power: [], aeroPower: [], thrust: [], torque: [], omega: [], pitchDeg: [], cp: [], ct: [], tsr: [], region: [],
    ratedWind: NaN, cpMax, tsrOpt, pitchFineDeg: fine };
  // rated wind speed: electrical power at fine pitch on the speed schedule equals rated power
  const g = (U) => eta * run(U, wAt(U), fine).power - Prated;
  {
    let lo = Math.max(cutIn, 1); let hi = Number.isFinite(cutOut) ? cutOut : 40;
    if (g(lo) < 0 && g(hi) > 0) {
      // the first crossing: march up in 1 m/s steps, then bisect
      let x = lo;
      while (x + 1 < hi && g(x + 1) < 0) x += 1;
      lo = x; hi = Math.min(hi, x + 1);
      for (let i = 0; i < 40; i++) { const mid = 0.5 * (lo + hi); if (g(mid) < 0) lo = mid; else hi = mid; }
      out.ratedWind = 0.5 * (lo + hi);
    }
  }
  for (const U of speeds) {
    let w = 0; let p = fine; let s = null; let region = 'off';
    if (U >= cutIn - 1e-9 && U <= cutOut + 1e-9 && U > 0) {
      w = wAt(U);
      s = run(U, w, fine);
      const wOpt = tsrOpt * U / R;
      region = wOpt < wMin ? '1.5' : (wOpt > wRated ? '2.5' : '2');
      if (eta * s.power > Prated) {
        region = '3';
        let lo = fine; let hi = fine + 2;
        let sh = run(U, w, hi);
        while (eta * sh.power > Prated && hi < 90) { lo = hi; hi += 2; sh = run(U, w, hi); }
        for (let i = 0; i < 40 && hi - lo > 1e-5; i++) {
          const mid = 0.5 * (lo + hi);
          if (eta * run(U, w, mid).power > Prated) lo = mid; else hi = mid;
        }
        p = 0.5 * (lo + hi);
        s = run(U, w, p);
      }
    }
    out.U.push(U); out.omega.push(w); out.pitchDeg.push(p); out.region.push(region);
    out.aeroPower.push(s ? s.power : 0); out.power.push(s ? Math.max(0, eta * s.power) : 0);
    out.thrust.push(s ? s.thrust : 0); out.torque.push(s ? s.torque : 0);
    out.cp.push(s ? s.Cp : 0); out.ct.push(s ? s.Ct : 0); out.tsr.push(s ? s.tsr : 0);
  }
  return out;
}
