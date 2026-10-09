/**
 * Lifting and lowering of a body through the water column on a cable: static tension, axial resonance of the
 * cable–object system, dynamic amplification under crane-tip heave with linearised quadratic drag, slack/snap
 * checks, heave-compensator transmissibility, pendulum period, horizontal offset in current, landing velocity and
 * the DNV-RP-N103 Simplified Method for splash-zone crossing.
 *
 * Conventions: depths and lengths in m (positive), forces in N (tension positive), z not used; amplitudes are
 * single amplitudes.
 */

import { brent } from '../math/roots.js';

const G0 = 9.80665;

function num(v, name, { min = -Infinity, max = Infinity, strictMin = false } = {}) {
  if (typeof v !== 'number' || !Number.isFinite(v) || v < min || v > max || (strictMin && v <= min)) {
    throw new RangeError(`lowering: argument ${name} = ${v} is outside its valid range`);
  }
  return v;
}

/**
 * Static hook load versus lowered depth: F = M·g − ρ·V·g + w·L (object weight in water plus submerged cable weight
 * paid out below the surface). Reference: DNV-RP-N103 §5.2.1. Validity: vertical cable, fully submerged object,
 * no current.
 * @param {object} o `{ mass [kg], volume [m³], cableWeight [N/m submerged], depths: number[] [m cable length in
 *   water], rho? = 1025, g? }`
 * @returns {{ objectWeightInWater: number, rows: Array<{ depth: number, tension: number }> }} [N]; [m, N at the top]
 */
export function staticTension({ mass, volume, cableWeight, depths, rho = 1025, g = G0 }) {
  num(mass, 'mass', { min: 0, strictMin: true }); num(volume, 'volume', { min: 0 }); num(cableWeight, 'cableWeight', { min: 0 });
  if (!Array.isArray(depths)) throw new RangeError('lowering: argument depths must be an array');
  const W = (mass - rho * volume) * g;
  return { objectWeightInWater: W, rows: depths.map((d, i) => ({ depth: num(d, `depths[${i}]`, { min: 0 }), tension: W + cableWeight * d })) };
}

/**
 * Axial natural period of the cable–object system.
 *
 * Exact (elastic rod of distributed mass m·L, stiffness EA, fixed at the top, end mass M' = M + A33):
 *   (ωL/c)·tan(ωL/c) = m·L/M',  c = sqrt(EA/m)   (Timoshenko, "Vibration Problems in Engineering", rod with end mass)
 * Lumped approximation (DNV-RP-N103 §5.3.3): T0 = 2π·sqrt((M + A33 + θ·m·L)/K), θ = 1/3, 1/K = 1/kc + L/EA
 * (kc = stiffness of crane tip / compensator / soft sling in series; exact equation solved with kc = ∞ unless kc
 * is given, in which case only the lumped value is returned for the combined system).
 * @param {object} o `{ mass [kg], addedMass? = 0 [kg heave added mass A33], EA [N], length [m], cableMass [kg/m],
 *   kc? [N/m], theta? = 1/3 }`
 * @returns {{ periodExact: number, periodLumped: number, omegaExact: number, omegaLumped: number, stiffness: number }}
 *   [s, s, rad/s, rad/s, N/m] (periodExact = NaN when kc is given)
 */
export function axialResonance({ mass, addedMass = 0, EA, length, cableMass, kc, theta = 1 / 3 }) {
  num(mass, 'mass', { min: 0, strictMin: true }); num(addedMass, 'addedMass', { min: 0 }); num(EA, 'EA', { min: 0, strictMin: true });
  num(length, 'length', { min: 0, strictMin: true }); num(cableMass, 'cableMass', { min: 0 });
  const Mp = mass + addedMass;
  const K = 1 / ((kc !== undefined ? 1 / num(kc, 'kc', { min: 0, strictMin: true }) : 0) + length / EA);
  const wl = Math.sqrt(K / (Mp + theta * cableMass * length));
  let we = NaN;
  if (kc === undefined) {
    if (cableMass === 0) we = Math.sqrt(EA / (length * Mp));
    else {
      const ratio = cableMass * length / Mp;
      const x = brent((q) => q * Math.tan(q) - ratio, 1e-12, Math.PI / 2 - 1e-12, { tol: 1e-14 });
      we = x * Math.sqrt(EA / cableMass) / length;
    }
  }
  return { periodExact: 2 * Math.PI / we, periodLumped: 2 * Math.PI / wl, omegaExact: we, omegaLumped: wl, stiffness: K };
}

/**
 * Dynamic response of the lowered object to harmonic crane-tip heave, single degree of freedom with
 * equivalent-linearised quadratic drag, and the slack-wire check.
 *
 * M'·ẍ + c·ẋ + K·(x − η) = 0 with M' = M + A33 + θ·m·L, η = ηa·sin ωt; c = (8/3π)·½ρ·Cd·Ap·ω·|x̂| (equivalent
 * linearisation of |ẋ|ẋ for a sinusoid, iterated). Dynamic line force amplitude Fd = K·|x̂ − η̂| =
 * K·ηa·|M'ω² − iωc| / |K − M'ω² + iωc|. Slack when Fd ≥ static tension at the object.
 * Reference: DNV-RP-N103 §5.3.5–5.3.7 (dynamic force and amplification, linearised damping); first-principles SDOF.
 * Validity: vertical cable, object far below the wave zone (no wave excitation), single harmonic.
 * @param {object} o `{ mass, addedMass? = 0 [kg], EA [N], length [m], cableMass [kg/m], kc? [N/m], heaveAmp [m],
 *   period [s], Cd? = 0 [-], area? = 0 [m² projected], rho? = 1025, linearDamping? = 0 [N·s/m], staticTension [N at
 *   the object], theta? = 1/3 }`
 * @returns {{ motionAmp: number, motionRatio: number, dynamicForce: number, daf: number, minTension: number,
 *   maxTension: number, slack: boolean, naturalPeriod: number, damping: number, iterations: number }}
 *   [m, -, N, -, N, N, -, s, N·s/m, -] (daf = (Fstatic + Fd)/Fstatic)
 */
export function heaveResponse(o) {
  const { mass, addedMass = 0, EA, length, cableMass, kc, heaveAmp, period, Cd = 0, area = 0, rho = 1025, linearDamping = 0, staticTension: Fs, theta = 1 / 3 } = o ?? {};
  num(heaveAmp, 'heaveAmp', { min: 0 }); num(period, 'period', { min: 0, strictMin: true }); num(Fs, 'staticTension', { min: 0, strictMin: true });
  num(Cd, 'Cd', { min: 0 }); num(area, 'area', { min: 0 }); num(linearDamping, 'linearDamping', { min: 0 });
  const res = axialResonance({ mass, addedMass, EA, length, cableMass, kc, theta });
  const K = res.stiffness, Mp = mass + addedMass + theta * cableMass * length, om = 2 * Math.PI / period;
  const kq = (8 / (3 * Math.PI)) * 0.5 * rho * Cd * area;
  let x = heaveAmp, c = linearDamping, it = 0;
  for (; it < 200; it++) {
    c = linearDamping + kq * om * x;
    const xn = heaveAmp * K / Math.hypot(K - Mp * om * om, om * c);
    if (Math.abs(xn - x) <= 1e-12 * (heaveAmp + 1e-30)) { x = xn; break; }
    x = 0.5 * (x + xn);
  }
  c = linearDamping + kq * om * x;
  const Fd = K * heaveAmp * Math.hypot(Mp * om * om, om * c) / Math.hypot(K - Mp * om * om, om * c);
  return {
    motionAmp: x, motionRatio: heaveAmp > 0 ? x / heaveAmp : 1, dynamicForce: Fd, daf: (Fs + Fd) / Fs,
    minTension: Fs - Fd, maxTension: Fs + Fd, slack: Fd >= Fs, naturalPeriod: res.periodLumped, damping: c, iterations: it + 1,
  };
}

/**
 * Passive / active heave compensator as a spring–damper stage between crane tip and load: motion
 * transmissibility of a base-excited mass, |x/η| = sqrt((k² + (cω)²) / ((k − Mω²)² + (cω)²)), and efficiency
 * 1 − |x/η|. An active stage is represented by the residual factor `activeResidual` (fraction of crane-tip motion
 * left by the active control, supplied by the caller) multiplying the passive transmissibility.
 * Reference: DNV-RP-N103 §5.5 (heave compensation, simplified model); classical base-excitation transmissibility
 * (Den Hartog, "Mechanical Vibrations"). Validity: linear spring and damper, rigid lift wire below the compensator.
 * @param {object} o `{ mass [kg incl. added mass], k [N/m], c [N·s/m], period [s], activeResidual? = 1 [-] }`
 * @returns {{ transmissibility: number, efficiency: number, naturalPeriod: number, dampingRatio: number }}
 *   [-, -, s, -]
 */
export function heaveCompensator({ mass, k, c, period, activeResidual = 1 }) {
  num(mass, 'mass', { min: 0, strictMin: true }); num(k, 'k', { min: 0, strictMin: true }); num(c, 'c', { min: 0 });
  num(period, 'period', { min: 0, strictMin: true }); num(activeResidual, 'activeResidual', { min: 0, max: 1 });
  const om = 2 * Math.PI / period;
  const tr = activeResidual * Math.sqrt((k * k + (c * om) ** 2) / ((k - mass * om * om) ** 2 + (c * om) ** 2));
  return { transmissibility: tr, efficiency: 1 - tr, naturalPeriod: 2 * Math.PI * Math.sqrt(mass / k), dampingRatio: c / (2 * Math.sqrt(k * mass)) };
}

/**
 * Pendulum period of a load hanging in air: simple pendulum T = 2π·sqrt(L/g); with a finite amplitude θ0 the
 * second-order correction T·(1 + θ0²/16) is applied (series expansion of the exact elliptic-integral period).
 * @param {object} o `{ length [m hoist length crane tip → load centre of gravity], amplitudeDeg? = 0, g? }`
 * @returns {number} period [s]
 */
export function pendulumPeriod({ length, amplitudeDeg = 0, g = G0 }) {
  num(length, 'length', { min: 0, strictMin: true }); num(amplitudeDeg, 'amplitudeDeg', { min: 0, max: 60 });
  const th = amplitudeDeg * Math.PI / 180;
  return 2 * Math.PI * Math.sqrt(length / g) * (1 + th * th / 16);
}

/**
 * Horizontal offset of a lowered object in a steady current: cable and object drag against the restoring effect
 * of the tension. From the object upwards (s' = height above the object), horizontal force
 * Fh(s') = FD0 + ∫q ds', tension T(s') = W + w·s', slope dx/ds' = Fh/T (small-angle form of DNV-RP-N103 §5.2.2);
 * integrated numerically for an arbitrary current profile. For uniform current the closed form is
 * x(L) = (q/w)·L + ((FD0 − q·W/w)/w)·ln(1 + w·L/W).
 * Validity: offsets small compared with the length (slope ≪ 1), steady current, W > 0.
 * @param {object} o `{ length [m], objectWeight [N weight in water W], cableWeight [N/m submerged w], cableDiameter [m],
 *   cableCd? = 1.2, objectDragArea [m² Cd·A of the object], current: number | (depth [m]) => U [m/s], rho? = 1025,
 *   nSteps? = 400 }`
 * @returns {{ offset: number, objectDrag: number, topAngleDeg: number, bottomAngleDeg: number, offsetProfile:
 *   { depth: number[], x: number[] } }} offset of the object relative to the crane tip [m]; drag [N]; cable angles from
 *   vertical [deg]; profile of horizontal position relative to the object versus depth
 */
export function currentOffset({ length, objectWeight, cableWeight, cableDiameter, cableCd = 1.2, objectDragArea, current, rho = 1025, nSteps = 400 }) {
  num(length, 'length', { min: 0, strictMin: true }); num(objectWeight, 'objectWeight', { min: 0, strictMin: true });
  num(cableWeight, 'cableWeight', { min: 0 }); num(cableDiameter, 'cableDiameter', { min: 0 }); num(objectDragArea, 'objectDragArea', { min: 0 });
  const U = typeof current === 'function' ? current : () => num(current, 'current', { min: 0 });
  const FD0 = 0.5 * rho * objectDragArea * U(length) ** 2;
  const q = (sp) => 0.5 * rho * cableCd * cableDiameter * U(length - sp) ** 2;
  const h = length / nSteps;
  let Fh = FD0, x = 0;
  const depth = [length], xs = [0];
  let slope0 = FD0 / objectWeight, slope = slope0;
  for (let k = 0; k < nSteps; k++) {
    const s0 = k * h, s1 = s0 + h;
    const Fh1 = Fh + 0.5 * (q(s0) + q(s1)) * h;
    const sl1 = Fh1 / (objectWeight + cableWeight * s1);
    x += 0.5 * (slope + sl1) * h;
    Fh = Fh1; slope = sl1;
    depth.push(length - s1); xs.push(x);
  }
  return { offset: x, objectDrag: FD0, topAngleDeg: Math.atan(slope) * 180 / Math.PI, bottomAngleDeg: Math.atan(slope0) * 180 / Math.PI, offsetProfile: { depth, x: xs } };
}

/**
 * Landing velocity estimate at the seabed: lowering speed plus the vertical velocity amplitude of the object,
 * v = vc + ω·|x̂| with |x̂| from heaveResponse (or supplied directly). Reference: DNV-RP-N103 §6.2 (landing on
 * seabed: impact velocity from hook lowering velocity and crane-tip motion transferred to the object).
 * @param {object} o `{ loweringSpeed [m/s], objectMotionAmp [m], period [s] }`
 * @returns {{ velocity: number, motionVelocity: number }} [m/s]
 */
export function landingVelocity({ loweringSpeed, objectMotionAmp, period }) {
  num(loweringSpeed, 'loweringSpeed', { min: 0 }); num(objectMotionAmp, 'objectMotionAmp', { min: 0 }); num(period, 'period', { min: 0, strictMin: true });
  const mv = 2 * Math.PI / period * objectMotionAmp;
  return { velocity: loweringSpeed + mv, motionVelocity: mv };
}

/* ------------------------------------------------------------------------------------------------------- */
/* DNV-RP-N103 Simplified Method (splash zone)                                                               */

/**
 * Characteristic crane-tip vertical motion from vessel heave, roll and pitch, DNV-RP-N103 §4.3.3:
 * ηct = sqrt(η3² + (b·sin η4)² + (l·sin η5)²), vct = 2π·sqrt((η3/T3)² + (b·sin η4/T4)² + (l·sin η5/T5)²),
 * act = 4π²·sqrt((η3/T3²)² + (b·sin η4/T4²)² + (l·sin η5/T5²)²).
 * @param {object} o `{ heaveAmp [m], rollAmpDeg? = 0, pitchAmpDeg? = 0 [deg], b? = 0 [m horizontal distance crane tip →
 *   roll axis], l? = 0 [m crane tip → pitch axis], T3, T4?, T5? [s natural periods of heave, roll, pitch] }`
 * @returns {{ eta: number, velocity: number, acceleration: number }} [m, m/s, m/s²] single amplitudes
 */
export function craneTipMotion({ heaveAmp, rollAmpDeg = 0, pitchAmpDeg = 0, b = 0, l = 0, T3, T4, T5 }) {
  num(heaveAmp, 'heaveAmp', { min: 0 }); num(T3, 'T3', { min: 0, strictMin: true });
  const r = b * Math.sin(rollAmpDeg * Math.PI / 180), p = l * Math.sin(pitchAmpDeg * Math.PI / 180);
  const t4 = r !== 0 ? num(T4, 'T4', { min: 0, strictMin: true }) : 1, t5 = p !== 0 ? num(T5, 'T5', { min: 0, strictMin: true }) : 1;
  return {
    eta: Math.hypot(heaveAmp, r, p),
    velocity: 2 * Math.PI * Math.hypot(heaveAmp / T3, r / t4, p / t5),
    acceleration: 4 * Math.PI ** 2 * Math.hypot(heaveAmp / T3 ** 2, r / t4 ** 2, p / t5 ** 2),
  };
}

/**
 * Characteristic wave kinematics of the DNV-RP-N103 Simplified Method (§4.3.2–4.3.4):
 * ζa = 0.9·Hs (operations of up to 30 minutes; pass `amplitudeFactor = 1.0` for longer operations),
 * vw = ζa·(2π/Tz)·e^(−4π²d/(Tz²·g)), aw = ζa·(2π/Tz)²·e^(−4π²d/(Tz²·g)), d = distance from the still water level to the
 * centre of the item (positive down). Stated validity of the method: 8.9·sqrt(Hs/g) ≤ Tz ≤ 13 s.
 * @param {object} o `{ Hs [m], Tz [s], d [m], amplitudeFactor? = 0.9, g? }`
 * @returns {{ zetaA: number, velocity: number, acceleration: number, tzInRange: boolean }} [m, m/s, m/s², -]
 */
export function waveKinematicsN103({ Hs, Tz, d, amplitudeFactor = 0.9, g = G0 }) {
  num(Hs, 'Hs', { min: 0 }); num(Tz, 'Tz', { min: 0, strictMin: true }); num(d, 'd', { min: 0 });
  const za = amplitudeFactor * Hs, om = 2 * Math.PI / Tz, dec = Math.exp(-4 * Math.PI ** 2 * d / (Tz * Tz * g));
  return { zetaA: za, velocity: za * om * dec, acceleration: za * om * om * dec, tzInRange: Tz >= 8.9 * Math.sqrt(Hs / g) && Tz <= 13 };
}

/**
 * Splash-zone hydrodynamic force and acceptance checks by the DNV-RP-N103 Simplified Method (§4.3–4.4), for one
 * item at one submergence level:
 *   slamming impact velocity   vs = vc + sqrt(vct² + vw²)                       (§4.3.5)
 *   slamming force             Fslam = ½·ρ·Cs·As·vs²
 *   varying buoyancy           Fρ = ρ·δV·g, δV = Ãw·sqrt(ζa² + ηct²)             (§4.3.6)
 *   mass force                 FM = sqrt([(M + A33)·act]² + [(ρV + A33)·aw]²)    (§4.3.7)
 *   drag force                 FD = ½·ρ·CD·Ap·vr², vr = vc + sqrt(vct² + vw²)    (§4.3.8)
 *   hydrodynamic force         Fhyd = sqrt((FD + Fslam)² + (FM − Fρ)²)           (§4.3.9)
 *   slack-sling criterion      Fhyd ≤ 0.9·Fstatic-min                           (§4.4.3)
 *   capacity check             Ftotal = Fstatic-max + Fhyd, DAFconv = Ftotal/(M·g) (§4.4.2, 4.4.4)
 * Coefficients Cs (not less than 3.0 for smooth circular cylinders and 5.0 otherwise, §4.3.5), CD (not less than
 * 2.5 in oscillatory flow unless documented, §4.3.8), A33 and the areas are arguments; so are the crane-tip
 * amplitudes (craneTipMotion) and wave kinematics (waveKinematicsN103).
 * @param {object} o `{ mass [kg], volume [m³ submerged at this level], addedMass [kg A33], slamArea? = 0 [m²],
 *   Cs? = 5, dragArea? = 0 [m² Ap], CD? = 2.5, waterplaneArea? = 0 [m² Ãw], vc? = 0.5 [m/s hook lowering velocity],
 *   etaCt [m], vCt [m/s], aCt [m/s²], zetaA [m], vW [m/s], aW [m/s²], staticMin? [N], staticMax? [N], rho? = 1025, g? }`
 *   staticMin / staticMax default to M·g − ρ·V·g
 * @returns {{ slamVelocity: number, Fslam: number, Frho: number, FM: number, FD: number, Fhyd: number,
 *   staticMin: number, staticMax: number, slackUc: number, slackOk: boolean, Ftotal: number, dafConv: number }}
 *   velocity [m/s]; forces [N]; slackUc = Fhyd/(0.9·Fstatic-min) [-]
 */
export function splashZoneN103(o) {
  const {
    mass, volume, addedMass, slamArea = 0, Cs = 5, dragArea = 0, CD = 2.5, waterplaneArea = 0, vc = 0.5,
    etaCt, vCt, aCt, zetaA, vW, aW, rho = 1025, g = G0,
  } = o ?? {};
  num(mass, 'mass', { min: 0, strictMin: true }); num(volume, 'volume', { min: 0 }); num(addedMass, 'addedMass', { min: 0 });
  for (const [v, n] of [[slamArea, 'slamArea'], [Cs, 'Cs'], [dragArea, 'dragArea'], [CD, 'CD'], [waterplaneArea, 'waterplaneArea'], [vc, 'vc'], [etaCt, 'etaCt'], [vCt, 'vCt'], [aCt, 'aCt'], [zetaA, 'zetaA'], [vW, 'vW'], [aW, 'aW']]) num(v, n, { min: 0 });
  const vs = vc + Math.hypot(vCt, vW);
  const Fslam = 0.5 * rho * Cs * slamArea * vs * vs;
  const Frho = rho * waterplaneArea * Math.hypot(zetaA, etaCt) * g;
  const FM = Math.hypot((mass + addedMass) * aCt, (rho * volume + addedMass) * aW);
  const FD = 0.5 * rho * CD * dragArea * vs * vs;
  const Fhyd = Math.hypot(FD + Fslam, FM - Frho);
  const stat = mass * g - rho * volume * g;
  const staticMin = o.staticMin ?? stat, staticMax = o.staticMax ?? stat;
  const Ftotal = staticMax + Fhyd;
  return {
    slamVelocity: vs, Fslam, Frho, FM, FD, Fhyd, staticMin, staticMax,
    slackUc: staticMin > 0 ? Fhyd / (0.9 * staticMin) : Infinity, slackOk: Fhyd <= 0.9 * staticMin,
    Ftotal, dafConv: Ftotal / (mass * g),
  };
}

/**
 * Snap load after a slack sling, DNV-RP-N103 §4.7.2: Fsnap = vsnap·sqrt(K·(M + A33)), vsnap = vff + C·vr, with
 * C = 1 for vff < 0.2·vr, C = cos[π·(vff/vr − 0.2)] for 0.2·vr ≤ vff ≤ 0.7·vr, C = 0 for vff > 0.7·vr; vr = vc +
 * sqrt(vct² + vw²). The free-fall velocity vff and the hoisting-system stiffness K are arguments.
 * @param {object} o `{ mass [kg], addedMass? = 0 [kg], K [N/m], vff [m/s], vr [m/s] }`
 * @returns {{ Fsnap: number, vSnap: number, C: number }} [N, m/s, -]
 */
export function snapLoadN103({ mass, addedMass = 0, K, vff, vr }) {
  num(mass, 'mass', { min: 0, strictMin: true }); num(addedMass, 'addedMass', { min: 0 }); num(K, 'K', { min: 0, strictMin: true });
  num(vff, 'vff', { min: 0 }); num(vr, 'vr', { min: 0, strictMin: true });
  const ratio = vff / vr;
  const C = ratio < 0.2 ? 1 : ratio <= 0.7 ? Math.cos(Math.PI * (ratio - 0.2)) : 0;
  const vSnap = vff + C * vr;
  return { Fsnap: vSnap * Math.sqrt(K * (mass + addedMass)), vSnap, C };
}
