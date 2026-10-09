/**
 * Ancillary equipment for flexible risers, umbilicals and cables: bend stiffeners, bend restrictors, bellmouths,
 * buoyancy modules and clamps, mid-water arches, pull-in through I/J-tubes, piggyback/saddle spacing, impact
 * protection, hang-off clamps and termination landing loads.
 */

function pos(v, name) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}
function fin(v, name) {
  if (!(typeof v === 'number' && Number.isFinite(v))) throw new RangeError(`${name} must be a finite number`);
  return v;
}
function nonneg(v, name) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v >= 0)) throw new RangeError(`${name} must be a non-negative finite number`);
  return v;
}

/**
 * Bend stiffener: tapered polyurethane cone on a riser/cable loaded by a tip tension T at angle theta to the
 * root axis.
 *
 * Model. Planar elastica of the composite beam, arc length s from the root (clamped, phi = 0):
 *     phi' = M / EI(s),   M' = -T sin(theta - phi)            (large angle; linearised: sin x -> x)
 *     EI(s) = EI_pipe + E_pu pi (d(s)^4 - d_i^4)/64,  d(s) = D_root + (D_tip - D_root) s/L  (linear taper)
 * Beyond the stiffener tip the pipe continues with EI_pipe to a straight line along the tension. The first
 * integral of the elastica there, EI phi'^2 / 2 = T (1 - cos(theta - phi)), gives the exact tip condition
 *     M(L) = 2 sqrt(T EI_pipe) sin((theta - phi(L))/2)        (linearised: sqrt(T EI_pipe)(theta - phi(L))).
 * The two-point problem is solved by shooting on the root curvature with RK4 and bisection. For constant EI and
 * small angles the solution is theta - phi = theta exp(-s sqrt(T/EI)), root curvature theta sqrt(T/EI).
 * Outputs follow from equilibrium: root shear T sin(theta), root axial force T cos(theta), root moment M(0)
 * (shared between pipe and stiffener in proportion to EI). The polyurethane outer-fibre strain is kappa d/2.
 * Reference: Boef & Out (1990) "Analysis of a flexible riser top connection with bend restrictor", OTC 6436
 * (beam formulation); elastica first integral (Love, Mathematical Theory of Elasticity, Ch. XIX).
 * Validity: linear elastic materials (PU modulus is strain-, rate- and temperature-dependent: use the secant
 * value for the load case), no gap between pipe and stiffener, pipe bending stiffness independent of curvature,
 * quasi-static, long pipe beyond the tip.
 * @param {object} o
 * @param {number} o.tension T, N (> 0)
 * @param {number} o.angle theta, rad (0 < theta < pi)
 * @param {number} o.length stiffener length L, m
 * @param {number} o.rootDiameter D_root, m
 * @param {number} o.tipDiameter D_tip, m (>= innerDiameter)
 * @param {number} o.innerDiameter d_i (pipe outside diameter), m
 * @param {number} o.Epu polyurethane modulus, Pa
 * @param {number} o.EIpipe pipe bending stiffness, N m^2
 * @param {boolean} [o.linear=false] use the linearised equations
 * @param {number} [o.nSteps=400] RK4 steps along the stiffener
 * @param {number} [o.mbr] pipe minimum bend radius, m (for the utilisation)
 * @param {number} [o.strainLimit] allowable PU strain, - (typical design values of a few per cent; argument)
 * @returns {{s:number[], curvature:number[], angle:number[], moment:number[], rootCurvature:number,
 *   tipCurvaturePipe:number, maxCurvature:number, maxCurvatureAt:number, rootMoment:number,
 *   rootMomentStiffener:number, rootShear:number, rootAxial:number, maxPuStrain:number,
 *   curvaturePerAngle:number, mbrUtilisation:number|null, strainUtilisation:number|null}} m, 1/m, rad, N m
 *   (total), 1/m, 1/m (in the bare pipe just beyond the tip), 1/m, m, N m, N m, N, N, -, 1/m per rad (maximum
 *   curvature / angle: transfer function for fatigue analysis in the small-angle range)
 */
export function bendStiffener({ tension, angle, length, rootDiameter, tipDiameter, innerDiameter, Epu, EIpipe, linear = false, nSteps = 400, mbr, strainLimit }) {
  const T = pos(tension, 'tension'), L = pos(length, 'length');
  if (!(angle > 0 && angle < Math.PI)) throw new RangeError('angle must be in (0, pi)');
  pos(innerDiameter, 'innerDiameter'); nonneg(Epu, 'Epu'); pos(EIpipe, 'EIpipe');
  if (!(tipDiameter >= innerDiameter)) throw new RangeError('tipDiameter must be at least innerDiameter');
  if (!(rootDiameter >= tipDiameter)) throw new RangeError('rootDiameter must be at least tipDiameter');
  if (!(Number.isInteger(nSteps) && nSteps >= 10)) throw new RangeError('nSteps must be an integer >= 10');
  const di4 = innerDiameter ** 4;
  const dia = (s) => rootDiameter + (tipDiameter - rootDiameter) * s / L;
  const EI = (s) => EIpipe + Epu * Math.PI * (dia(s) ** 4 - di4) / 64;
  const load = linear ? (x) => x : Math.sin;
  const h = L / nSteps;
  const phi = new Float64Array(nSteps + 1), M = new Float64Array(nSteps + 1);
  const march = (k0) => {
    let p = 0, m = EI(0) * k0;
    phi[0] = 0; M[0] = m;
    for (let i = 0; i < nSteps; i++) {
      const s = i * h, eiA = EI(s), eiB = EI(s + h / 2), eiC = EI(s + h);
      const k1p = m / eiA, k1m = -T * load(angle - p);
      const k2p = (m + h / 2 * k1m) / eiB, k2m = -T * load(angle - (p + h / 2 * k1p));
      const k3p = (m + h / 2 * k2m) / eiB, k3m = -T * load(angle - (p + h / 2 * k2p));
      const k4p = (m + h * k3m) / eiC, k4m = -T * load(angle - (p + h * k3p));
      p += h / 6 * (k1p + 2 * k2p + 2 * k3p + k4p);
      m += h / 6 * (k1m + 2 * k2m + 2 * k3m + k4m);
      phi[i + 1] = p; M[i + 1] = m;
    }
    const d = angle - p;
    return m - (linear ? Math.sqrt(T * EIpipe) * d : 2 * Math.sqrt(T * EIpipe) * Math.sin(d / 2));
  };
  let lo = 0, hi = angle * Math.sqrt(T / EIpipe);
  for (let i = 0; i < 80 && !(march(hi) > 0); i++) hi *= 2;
  if (!(march(hi) > 0)) throw new Error('Bend stiffener shooting failed to bracket the root curvature');
  for (let i = 0; i < 200; i++) {
    const mid = 0.5 * (lo + hi);
    if (march(mid) > 0) hi = mid; else lo = mid;
    if (hi - lo <= 1e-15 * hi) break;
  }
  march(0.5 * (lo + hi));
  const s = [], curvature = [], ang = [], mom = [];
  let kMax = 0, sMax = 0, eMax = 0;
  for (let i = 0; i <= nSteps; i++) {
    const x = i * h, k = M[i] / EI(x);
    s.push(x); curvature.push(k); ang.push(phi[i]); mom.push(M[i]);
    if (Math.abs(k) > kMax) { kMax = Math.abs(k); sMax = x; }
    eMax = Math.max(eMax, Math.abs(k) * dia(x) / 2);
  }
  const tipPipe = M[nSteps] / EIpipe;
  if (Math.abs(tipPipe) > kMax) { kMax = Math.abs(tipPipe); sMax = L; }
  return {
    s, curvature, angle: ang, moment: mom, rootCurvature: curvature[0], tipCurvaturePipe: tipPipe, maxCurvature: kMax, maxCurvatureAt: sMax,
    rootMoment: M[0], rootMomentStiffener: M[0] * (1 - EIpipe / EI(0)), rootShear: T * Math.sin(angle), rootAxial: T * Math.cos(angle),
    maxPuStrain: eMax, curvaturePerAngle: kMax / angle,
    mbrUtilisation: mbr === undefined ? null : kMax * pos(mbr, 'mbr'),
    strainUtilisation: strainLimit === undefined ? null : eMax / pos(strainLimit, 'strainLimit'),
  };
}

/**
 * Bend stiffener sizing loop: scans the supplied lengths and a grid of root diameters and returns the feasible
 * design of least polyurethane volume for which the maximum pipe curvature does not exceed 1/mbr and the PU
 * strain does not exceed its limit, for every load case (tension, angle).
 * PU volume of the cone frustum with a cylindrical bore: V = pi L (D_r^2 + D_r D_t + D_t^2)/12 - pi d_i^2 L/4.
 * Fatigue hook: run {@link bendStiffener} on the selected design; its `curvaturePerAngle` is the transfer
 * function from end rotation to pipe curvature for the fatigue load bins.
 * @param {object} o
 * @param {{tension:number, angle:number}[]} o.loadCases N, rad
 * @param {number[]} o.lengths candidate lengths, m
 * @param {number} o.rootDiameterMax m
 * @param {number} o.tipDiameter m
 * @param {number} o.innerDiameter m
 * @param {number} o.Epu Pa
 * @param {number} o.EIpipe N m^2
 * @param {number} o.mbr m
 * @param {number} [o.strainLimit] -
 * @param {number} [o.nDiameters=25] grid points between tipDiameter and rootDiameterMax
 * @returns {{feasible:boolean, length:number|null, rootDiameter:number|null, volume:number|null,
 *   maxCurvature:number|null, maxPuStrain:number|null, utilisation:number|null}} m, m, m^3, 1/m, -, -
 */
export function sizeBendStiffener({ loadCases, lengths, rootDiameterMax, tipDiameter, innerDiameter, Epu, EIpipe, mbr, strainLimit, nDiameters = 25 }) {
  if (!Array.isArray(loadCases) || loadCases.length === 0) throw new RangeError('loadCases must be a non-empty array');
  if (!Array.isArray(lengths) || lengths.length === 0) throw new RangeError('lengths must be a non-empty array');
  pos(mbr, 'mbr');
  if (!(rootDiameterMax > tipDiameter)) throw new RangeError('rootDiameterMax must exceed tipDiameter');
  let best = null;
  for (const L of lengths) {
    for (let j = 0; j < nDiameters; j++) {
      const Dr = tipDiameter + (rootDiameterMax - tipDiameter) * (j + 1) / nDiameters;
      let k = 0, e = 0;
      for (const lc of loadCases) {
        const r = bendStiffener({ tension: lc.tension, angle: lc.angle, length: L, rootDiameter: Dr, tipDiameter, innerDiameter, Epu, EIpipe, nSteps: 200 });
        k = Math.max(k, r.maxCurvature); e = Math.max(e, r.maxPuStrain);
      }
      const ut = Math.max(k * mbr, strainLimit === undefined ? 0 : e / strainLimit);
      if (ut <= 1) {
        const volume = Math.PI * L * (Dr * Dr + Dr * tipDiameter + tipDiameter * tipDiameter) / 12 - Math.PI * innerDiameter ** 2 * L / 4;
        if (!best || volume < best.volume) best = { feasible: true, length: L, rootDiameter: Dr, volume, maxCurvature: k, maxPuStrain: e, utilisation: ut };
        break; // larger roots at this length only add volume
      }
    }
  }
  return best ?? { feasible: false, length: null, rootDiameter: null, volume: null, maxCurvature: null, maxPuStrain: null, utilisation: null };
}

/**
 * Bend restrictor (vertebrae) locking geometry and moment check.
 * Elements of pitch p that can each rotate by delta before locking form a polygon whose joints lie on a circle
 * of radius R_lock = p / (2 sin(delta/2)) (exact chord geometry; ~ p/delta). n elements accommodate a total
 * angle n delta. A locked restrictor carrying a tension T that leaves at angle theta sees at its root the moment
 * of T about the root, M = T R_lock (1 - cos(theta)) (perpendicular distance from the root to the tangent at the
 * end of a circular arc of angle theta), and the shear T sin(theta).
 * Reference: geometry and statics; API RP 17B §4 / API Spec 17L1 (bend restrictors: locking radius >= MBR,
 * structural capacity >= design moment).
 * Validity: theta <= n delta (otherwise the pipe over-bends beyond the restrictor: flagged).
 * @param {object} o
 * @param {number} o.pitch element pitch p, m
 * @param {number} o.lockAngle delta per element, rad
 * @param {number} o.nElements number of elements
 * @param {number} [o.mbr] pipe minimum bend radius, m
 * @param {number} [o.tension] N
 * @param {number} [o.angle] exit angle theta, rad
 * @param {number} [o.momentCapacity] N m
 * @param {number} [o.shearCapacity] N
 * @returns {{lockingRadius:number, totalAngle:number, lockedLength:number, mbrUtilisation:number|null,
 *   coversAngle:boolean|null, rootMoment:number|null, rootShear:number|null, momentUtilisation:number|null,
 *   shearUtilisation:number|null}} m, rad, m, -, -, N m, N, -, -
 */
export function bendRestrictor({ pitch, lockAngle, nElements, mbr, tension, angle, momentCapacity, shearCapacity }) {
  pos(pitch, 'pitch'); pos(nElements, 'nElements');
  if (!(lockAngle > 0 && lockAngle < Math.PI)) throw new RangeError('lockAngle must be in (0, pi)');
  const R = pitch / (2 * Math.sin(lockAngle / 2));
  const out = { lockingRadius: R, totalAngle: nElements * lockAngle, lockedLength: nElements * pitch, mbrUtilisation: mbr === undefined ? null : pos(mbr, 'mbr') / R, coversAngle: null, rootMoment: null, rootShear: null, momentUtilisation: null, shearUtilisation: null };
  if (tension !== undefined && angle !== undefined) {
    nonneg(tension, 'tension'); nonneg(angle, 'angle');
    const th = Math.min(angle, out.totalAngle);
    out.coversAngle = angle <= out.totalAngle * (1 + 1e-12);
    out.rootMoment = tension * R * (1 - Math.cos(th));
    out.rootShear = tension * Math.sin(angle);
    if (momentCapacity !== undefined) out.momentUtilisation = out.rootMoment / pos(momentCapacity, 'momentCapacity');
    if (shearCapacity !== undefined) out.shearUtilisation = out.rootShear / pos(shearCapacity, 'shearCapacity');
  }
  return out;
}

/**
 * Bellmouth with a constant-curvature (circular) profile of radius R for a maximum exit angle theta:
 * axial length R sin(theta), radial flare R (1 - cos(theta)) per side, exit diameter = throat diameter +
 * 2 R (1 - cos(theta)), contact line load on the profile T / R.
 * Reference: geometry; capstan statics for the line load. R should not be less than the product MBR.
 * @param {object} o
 * @param {number} o.radius profile radius R, m
 * @param {number} o.maxAngle theta, rad (0, pi/2]
 * @param {number} o.throatDiameter m
 * @param {number} [o.mbr] m
 * @param {number} [o.tension] N
 * @returns {{length:number, flare:number, exitDiameter:number, arcLength:number, mbrUtilisation:number|null,
 *   lineLoad:number|null}} m, m, m, m, -, N/m
 */
export function bellmouth({ radius, maxAngle, throatDiameter, mbr, tension }) {
  pos(radius, 'radius'); pos(throatDiameter, 'throatDiameter');
  if (!(maxAngle > 0 && maxAngle <= Math.PI / 2)) throw new RangeError('maxAngle must be in (0, pi/2]');
  const flare = radius * (1 - Math.cos(maxAngle));
  return {
    length: radius * Math.sin(maxAngle), flare, exitDiameter: throatDiameter + 2 * flare, arcLength: radius * maxAngle,
    mbrUtilisation: mbr === undefined ? null : pos(mbr, 'mbr') / radius,
    lineLoad: tension === undefined ? null : nonneg(tension, 'tension') / radius,
  };
}

/**
 * GENERIC, editable table of syntactic-foam density against depth rating (kg/m^3): indicative round values
 * showing the trend that deeper-rated foams are denser; NOT supplier data.
 */
export const GENERIC_FOAM_DENSITY = Object.freeze({
  generic: true,
  depth: Object.freeze([500, 1000, 2000, 3000, 4000]),
  density: Object.freeze([350, 400, 480, 550, 620]),
});

/**
 * Distributed buoyancy modules: net buoyancy per module, number and pitch for a required distributed buoyancy,
 * and clamp slip check.
 *   foam volume V = pi (Do^2 - Di^2) L / 4;  net buoyancy B = [(rho_w - rho_foam) V] g (1 - lossAllowance) - W_hw
 *   (W_hw = submerged weight of clamp and hardware);  pitch = B / b_req;  count = ceil(b_req L_section / B).
 * Foam density is given directly or interpolated linearly in the (editable) depth-rating table.
 * The water-ingress / hydrostatic-compression buoyancy loss over the service life is an argument (suppliers
 * quote a guaranteed end-of-life loss).
 * Clamp slip: the clamp must hold the module's force component along the riser: capacity mu N_clamp against
 * demand = axialLoad (argument; by default the full net buoyancy, i.e. a vertical riser) x safety factor.
 * Reference: Archimedes; API RP 17B §4 / API Spec 17L1 (buoyancy modules and clamps) in intent.
 * @param {object} o
 * @param {number} o.outerDiameter Do, m
 * @param {number} o.innerDiameter Di (bore for the riser and clamp recess), m
 * @param {number} o.moduleLength L, m
 * @param {number} [o.foamDensity] kg/m^3
 * @param {number} [o.depthRating] m (used with the table when foamDensity is omitted)
 * @param {object} [o.table=GENERIC_FOAM_DENSITY] { depth: [], density: [] }
 * @param {number} [o.hardwareWeight=0] submerged weight of clamp/hardware per module, N
 * @param {number} [o.lossAllowance=0] fractional buoyancy loss, [0, 1)
 * @param {number} [o.requiredBuoyancy] distributed net buoyancy b_req, N/m
 * @param {number} [o.sectionLength] buoyant section length, m
 * @param {number} [o.clampForce] total radial clamping force N_clamp, N
 * @param {number} [o.mu] clamp/riser friction coefficient
 * @param {number} [o.axialLoad] slip demand along the riser, N (default: net buoyancy)
 * @param {number} [o.slipSafetyFactor=1.5]
 * @param {number} [o.rhoW=1025]
 * @param {number} [o.g=9.80665]
 * @returns {{foamDensity:number, foamVolume:number, netBuoyancy:number, pitch:number|null, count:number|null,
 *   providedBuoyancy:number|null, slipCapacity:number|null, slipUtilisation:number|null,
 *   genericDensity:boolean}} kg/m^3, m^3, N, m, -, N/m, N, -
 */
export function buoyancyModules({ outerDiameter, innerDiameter, moduleLength, foamDensity, depthRating, table = GENERIC_FOAM_DENSITY, hardwareWeight = 0, lossAllowance = 0, requiredBuoyancy, sectionLength, clampForce, mu, axialLoad, slipSafetyFactor = 1.5, rhoW = 1025, g = 9.80665 }) {
  pos(outerDiameter, 'outerDiameter'); pos(innerDiameter, 'innerDiameter'); pos(moduleLength, 'moduleLength');
  if (!(outerDiameter > innerDiameter)) throw new RangeError('outerDiameter must exceed innerDiameter');
  nonneg(hardwareWeight, 'hardwareWeight'); pos(rhoW, 'rhoW'); pos(g, 'g');
  if (!(lossAllowance >= 0 && lossAllowance < 1)) throw new RangeError('lossAllowance must be in [0, 1)');
  let rho = foamDensity, genericDensity = false;
  if (rho === undefined) {
    pos(depthRating, 'depthRating');
    const xs = table.depth, ys = table.density;
    if (!Array.isArray(xs) && !ArrayBuffer.isView(xs)) throw new RangeError('table.depth must be an array');
    if (depthRating > xs[xs.length - 1]) throw new RangeError('depthRating exceeds the foam table');
    let i = 0;
    while (i < xs.length - 2 && depthRating > xs[i + 1]) i++;
    rho = depthRating <= xs[0] ? ys[0] : ys[i] + (ys[i + 1] - ys[i]) * (depthRating - xs[i]) / (xs[i + 1] - xs[i]);
    genericDensity = table.generic === true;
  }
  pos(rho, 'foamDensity');
  const V = Math.PI * (outerDiameter ** 2 - innerDiameter ** 2) / 4 * moduleLength;
  const B = (rhoW - rho) * V * g * (1 - lossAllowance) - hardwareWeight;
  const out = { foamDensity: rho, foamVolume: V, netBuoyancy: B, pitch: null, count: null, providedBuoyancy: null, slipCapacity: null, slipUtilisation: null, genericDensity };
  if (requiredBuoyancy !== undefined) {
    pos(requiredBuoyancy, 'requiredBuoyancy');
    if (!(B > 0)) throw new RangeError('module has no net buoyancy');
    out.pitch = B / requiredBuoyancy;
    if (sectionLength !== undefined) {
      out.count = Math.ceil(requiredBuoyancy * pos(sectionLength, 'sectionLength') / B - 1e-9);
      out.providedBuoyancy = out.count * B / sectionLength;
    }
  }
  if (clampForce !== undefined && mu !== undefined) {
    out.slipCapacity = pos(mu, 'mu') * pos(clampForce, 'clampForce');
    out.slipUtilisation = pos(slipSafetyFactor, 'slipSafetyFactor') * Math.abs(axialLoad ?? B) / out.slipCapacity;
  }
  return out;
}

/**
 * Mid-water arch static balance. The arch (net buoyancy B, upward) carries two riser legs whose tensions T1 and
 * T2 pull away from the arch at angles a1, a2 below the horizontal on opposite sides (leg 1 towards -x, leg 2
 * towards +x). The tether supplies the remaining force:
 *     tether vertical (downward on the arch)   Tv = B - T1 sin(a1) - T2 sin(a2)
 *     tether horizontal                        Th = T1 cos(a1) - T2 cos(a2)   (towards +x on the arch)
 *     tether tension = sqrt(Tv^2 + Th^2), inclination from vertical = atan2(|Th|, Tv).
 * Tv <= 0 means the tether is slack (insufficient buoyancy).
 * Reference: statics of a point (force polygon).
 * @param {object} o
 * @param {number} o.netBuoyancy B, N
 * @param {number} o.T1 N
 * @param {number} o.angle1 a1, rad below horizontal
 * @param {number} o.T2 N
 * @param {number} o.angle2 a2, rad below horizontal
 * @param {number} [o.tetherCapacity] N
 * @param {number} [o.nTethers=1] tethers sharing the load equally
 * @returns {{tetherVertical:number, tetherHorizontal:number, tetherTension:number, tetherAngle:number,
 *   perTether:number, slack:boolean, utilisation:number|null}} N, N, N, rad from vertical, N, -, -
 */
export function midWaterArch({ netBuoyancy, T1, angle1, T2, angle2, tetherCapacity, nTethers = 1 }) {
  fin(netBuoyancy, 'netBuoyancy'); nonneg(T1, 'T1'); nonneg(T2, 'T2'); fin(angle1, 'angle1'); fin(angle2, 'angle2'); pos(nTethers, 'nTethers');
  const Tv = netBuoyancy - T1 * Math.sin(angle1) - T2 * Math.sin(angle2);
  const Th = T1 * Math.cos(angle1) - T2 * Math.cos(angle2);
  const Tt = Math.hypot(Tv, Th);
  return {
    tetherVertical: Tv, tetherHorizontal: Th, tetherTension: Tt, tetherAngle: Math.atan2(Math.abs(Th), Tv), perTether: Tt / nTethers, slack: Tv <= 0,
    utilisation: tetherCapacity === undefined ? null : Tt / nTethers / pos(tetherCapacity, 'tetherCapacity'),
  };
}

/**
 * Capstan (Euler-Eytelwein) relation for a line wrapped over a curved surface with friction:
 * T_tight = T_slack exp(mu theta).
 * Reference: Euler (1769) / Eytelwein (1808).
 * @param {object} o
 * @param {number} o.tensionIn slack-side tension, N
 * @param {number} o.mu friction coefficient
 * @param {number} o.wrapAngle theta, rad
 * @returns {{tensionOut:number, ratio:number}} N, -
 */
export function capstan({ tensionIn, mu, wrapAngle }) {
  nonneg(tensionIn, 'tensionIn'); nonneg(mu, 'mu'); nonneg(wrapAngle, 'wrapAngle');
  const ratio = Math.exp(mu * wrapAngle);
  return { tensionOut: tensionIn * ratio, ratio };
}

/**
 * Pull-in tension through an I-tube / J-tube or along a route made of straight sections and bends, marching
 * from the tail (back tension) to the pulling head:
 *   straight of length L at inclination beta above horizontal (pulling upward positive), submerged weight w,
 *     friction mu:  T2 = T1 + w L (sin(beta) + mu cos(beta))
 *   bend of angle theta with friction mu (weight within the bend neglected):  T2 = T1 exp(mu theta)
 *     and the sidewall line load there is T2 / R when the bend radius R is given.
 * Reference: standard cable-pulling equations (capstan + inclined-plane friction), e.g. Rifenburg (1953) AIEE
 * Trans. 72, "Pipe-line design for pipe-type feeders".
 * Validity: no bending-stiffness (bend entry) forces - significant for stiff products in tight J-tube bends;
 * no dynamic amplification.
 * @param {object} o
 * @param {number} [o.backTension=0] N
 * @param {({type:'straight', length:number, inclination?:number, mu:number, weight:number}|
 *   {type:'bend', angle:number, mu:number, radius?:number})[]} o.segments in pulling order; length m,
 *   inclination rad, weight N/m, angle rad, radius m
 * @param {number} [o.capacity] allowable pulling tension, N
 * @returns {{pullTension:number, tensions:number[], maxSidewall:number, utilisation:number|null}} N, N after
 *   each segment, N/m, -
 */
export function pullInTension({ backTension = 0, segments, capacity }) {
  nonneg(backTension, 'backTension');
  if (!Array.isArray(segments)) throw new RangeError('segments must be an array');
  let T = backTension, sw = 0;
  const tensions = segments.map((g, i) => {
    nonneg(g.mu, `segments[${i}].mu`);
    if (g.type === 'straight') {
      const b = fin(g.inclination ?? 0, `segments[${i}].inclination`);
      T = Math.max(T + nonneg(g.weight, `segments[${i}].weight`) * nonneg(g.length, `segments[${i}].length`) * (Math.sin(b) + g.mu * Math.cos(b)), 0);
    } else if (g.type === 'bend') {
      T = capstan({ tensionIn: T, mu: g.mu, wrapAngle: g.angle }).tensionOut;
      if (g.radius !== undefined) sw = Math.max(sw, T / pos(g.radius, `segments[${i}].radius`));
    } else throw new RangeError(`segments[${i}].type must be 'straight' or 'bend'`);
    return T;
  });
  return { pullTension: T, tensions, maxSidewall: sw, utilisation: capacity === undefined ? null : T / pos(capacity, 'capacity') };
}

/**
 * Spacing of clamps / piggyback saddles along a supported line (the piggybacked product spans between saddles as
 * a continuous beam under a uniform transverse load q, e.g. weight plus hydrodynamic load):
 *   bending:   M = q L^2 / 12 at the supports of a multi-span beam  ->  L <= sqrt(12 M_allow / q)
 *   clamp load: reaction q L per clamp                              ->  L <= P_clamp / q
 *   deflection (optional): 5 q L^4/(384 EI) for a single span is conservative for an interior span
 *                                                                   ->  L <= (384 EI delta/(5 q))^(1/4)
 * Reference: beam formulas (Roark, Table 8.1). Vortex-induced vibration of the span is a separate check.
 * @param {object} o
 * @param {number} o.load q, N/m
 * @param {number} [o.momentAllowable] N m
 * @param {number} [o.clampCapacity] N
 * @param {number} [o.EI] N m^2
 * @param {number} [o.deflectionLimit] m
 * @param {number} [o.length] line length for the clamp count, m
 * @returns {{spacing:number, governing:string, spacingMoment:number|null, spacingClamp:number|null,
 *   spacingDeflection:number|null, count:number|null}} m
 */
export function clampSpacing({ load, momentAllowable, clampCapacity, EI, deflectionLimit, length }) {
  const q = pos(load, 'load');
  const sM = momentAllowable === undefined ? null : Math.sqrt(12 * pos(momentAllowable, 'momentAllowable') / q);
  const sC = clampCapacity === undefined ? null : pos(clampCapacity, 'clampCapacity') / q;
  const sD = EI === undefined || deflectionLimit === undefined ? null : (384 * pos(EI, 'EI') * pos(deflectionLimit, 'deflectionLimit') / (5 * q)) ** 0.25;
  const c = [['bending', sM], ['clamp capacity', sC], ['deflection', sD]].filter((x) => x[1] !== null).sort((a, b) => a[1] - b[1]);
  if (c.length === 0) throw new RangeError('give at least one of momentAllowable, clampCapacity or EI with deflectionLimit');
  return { spacing: c[0][1], governing: c[0][0], spacingMoment: sM, spacingClamp: sC, spacingDeflection: sD, count: length === undefined ? null : Math.ceil(pos(length, 'length') / c[0][1] - 1e-9) + 1 };
}

/**
 * Impact / abrasion protection sleeves (polyurethane half-shells such as Uraduct-type products): energy check
 * against the rated absorption for the selected wall thickness (supplier test value, argument) and take-off.
 * Dropped-object energy if not given: E = m g h in air, or 0.5 m v^2 for a terminal velocity in water.
 * Reference: energy balance; DNV-RP-F107 impact-energy approach in intent.
 * @param {object} o
 * @param {number} [o.impactEnergy] J
 * @param {number} [o.mass] kg (with dropHeight or velocity)
 * @param {number} [o.dropHeight] m
 * @param {number} [o.velocity] m/s
 * @param {number} o.ratedEnergy energy absorbed without damage to the product, J
 * @param {number} [o.length] protected length, m
 * @param {number} [o.shellLength] m
 * @param {number} [o.g=9.80665]
 * @returns {{impactEnergy:number, utilisation:number, shells:number|null}} J, -, number of half-shell pairs
 */
export function impactProtection({ impactEnergy, mass, dropHeight, velocity, ratedEnergy, length, shellLength, g = 9.80665 }) {
  pos(ratedEnergy, 'ratedEnergy');
  let E = impactEnergy;
  if (E === undefined) {
    pos(mass, 'mass');
    E = velocity !== undefined ? 0.5 * mass * nonneg(velocity, 'velocity') ** 2 : mass * g * nonneg(dropHeight, 'dropHeight');
  }
  nonneg(E, 'impactEnergy');
  return { impactEnergy: E, utilisation: E / ratedEnergy, shells: length === undefined ? null : Math.ceil(pos(length, 'length') / pos(shellLength, 'shellLength') - 1e-9) };
}

/**
 * Hang-off friction clamp check: bolted split clamp gripping the armour/sheath over a length L.
 *   total clamping force N = nBolts x boltPreload (each bolt squeezes the product across the split)
 *   slip capacity = mu x N x faces / safetyFactor against the hang-off tension. The total radial normal force
 *   on the product lies between 2 N (two opposite line contacts, faces = 2, default and conservative) and pi N
 *   (uniform contact pressure round the circumference, faces = pi)
 *   mean radial contact pressure p = N / (D L) (bolt force over the projected area), against the allowable
 *   crush pressure of the product (argument).
 * Reference: friction-grip statics. Preload scatter, creep relaxation of polymer liners and thermal effects are
 * covered by the safety factor (argument).
 * @param {object} o
 * @param {number} o.tension hang-off tension, N
 * @param {number} o.nBolts
 * @param {number} o.boltPreload N per bolt
 * @param {number} o.mu friction coefficient
 * @param {number} o.diameter product diameter, m
 * @param {number} o.clampLength m
 * @param {number} [o.faces=2]
 * @param {number} [o.safetyFactor=2]
 * @param {number} [o.allowablePressure] Pa
 * @returns {{clampingForce:number, slipCapacity:number, slipUtilisation:number, contactPressure:number,
 *   pressureUtilisation:number|null}} N, N (factored), -, Pa, -
 */
export function hangOffClamp({ tension, nBolts, boltPreload, mu, diameter, clampLength, faces = 2, safetyFactor = 2, allowablePressure }) {
  nonneg(tension, 'tension'); pos(nBolts, 'nBolts'); pos(boltPreload, 'boltPreload'); pos(mu, 'mu'); pos(diameter, 'diameter'); pos(clampLength, 'clampLength');
  pos(faces, 'faces'); pos(safetyFactor, 'safetyFactor');
  const N = nBolts * boltPreload;
  const cap = mu * N * faces / safetyFactor;
  const p = N / (diameter * clampLength);
  return { clampingForce: N, slipCapacity: cap, slipUtilisation: tension / cap, contactPressure: p, pressureUtilisation: allowablePressure === undefined ? null : p / pos(allowablePressure, 'allowablePressure') };
}

/**
 * Landing load of a subsea termination (UTA / SUTU / cable end module) on the seabed or its foundation (hook):
 * single-degree-of-freedom impact of mass m (including added mass) with submerged weight W arriving at velocity
 * v on a linear contact stiffness k. Energy balance 0.5 m v^2 + W x = 0.5 k x^2 gives the peak reaction
 *     F = k x = W + sqrt(W^2 + k m v^2),
 * which tends to 2 W for v = 0 (suddenly applied load) and to v sqrt(k m) for hard impacts. `dynamicForce` is
 * F - W.
 * Reference: energy method for impact (landing on seabed, DNV-RP-N103 in intent).
 * Validity: linear stiffness, no soil plastic penetration or water cushioning (conservative for soft seabeds).
 * @param {object} o
 * @param {number} o.mass structure mass, kg
 * @param {number} [o.addedMass=0] kg
 * @param {number} o.submergedWeight N
 * @param {number} o.velocity landing velocity, m/s
 * @param {number} o.stiffness contact stiffness, N/m
 * @param {number} [o.capacity] allowable landing load, N
 * @returns {{dynamicForce:number, peakLoad:number, daf:number, kineticEnergy:number,
 *   utilisation:number|null}} N, N, -, J, -
 */
export function terminationLandingLoad({ mass, addedMass = 0, submergedWeight, velocity, stiffness, capacity }) {
  pos(mass, 'mass'); nonneg(addedMass, 'addedMass'); pos(submergedWeight, 'submergedWeight'); nonneg(velocity, 'velocity'); pos(stiffness, 'stiffness');
  const m = mass + addedMass;
  const peak = submergedWeight + Math.sqrt(submergedWeight ** 2 + stiffness * m * velocity ** 2);
  const F = peak - submergedWeight;
  return { dynamicForce: F, peakLoad: peak, daf: peak / submergedWeight, kineticEnergy: 0.5 * m * velocity ** 2, utilisation: capacity === undefined ? null : peak / pos(capacity, 'capacity') };
}
