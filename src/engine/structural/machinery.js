/**
 * Machinery vibration: rotor unbalance forcing, foundation response, vibration isolation and resonance checks.
 *
 *  - Residual unbalance and rotating force from the balance quality grade G = e·Ω (ISO 21940-11, formerly ISO 1940-1).
 *  - Harmonic response of a single-degree-of-freedom machine–foundation system, force transmissibility and isolation
 *    efficiency (Den Hartog, Mechanical Vibrations; Harris' Shock and Vibration Handbook ch. 2).
 *  - Two-degree-of-freedom model: machine on isolators on a foundation block (or deck) on its own support.
 *  - Stiffness and damping of a rigid circular footing on an elastic half-space (Richart, Hall & Woods, Vibrations of
 *    Soils and Foundations, 1970).
 *  - Separation margin between excitation and natural frequencies, vibration severity zones against supplied limits,
 *    conversion between harmonic displacement, velocity and acceleration, Campbell diagram data with crossings.
 *
 * SI units (m, kg, s, N, rad/s). Arguments whose name ends in `Rpm` are in revolutions per minute, `Hz` in hertz,
 * `MmS` in mm/s (the unit in which balance grades are designated).
 */
import { brent } from '../math/roots.js';

const G0 = 9.80665;

function pos(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}
function nonNeg(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v >= 0)) throw new RangeError(`${name} must be a non-negative finite number`);
  return v;
}

function speedOf({ omega, speedRpm, frequencyHz }) {
  if (omega !== undefined) return pos('omega', omega);
  if (speedRpm !== undefined) return pos('speedRpm', speedRpm) * Math.PI / 30;
  if (frequencyHz !== undefined) return pos('frequencyHz', frequencyHz) * 2 * Math.PI;
  throw new RangeError('one of omega, speedRpm or frequencyHz is required');
}

/* ------------------------------------------------------------------ unbalance */

/**
 * Rotating unbalance force of a rotor. With a balance quality grade G (ISO 21940-11: the product of the permissible
 * residual specific unbalance and the service speed, e_per·Ω = G, designated in mm/s, e.g. G 6.3, G 2.5) the
 * permissible eccentricity is e = G/Ω and the force amplitude F = m e Ω² = m G Ω. Alternatively give the eccentricity
 * directly. A service factor may be applied to allow for deterioration in operation (project or code value).
 * @param {object} o
 * @param {number} o.mass rotor mass [kg]
 * @param {number} [o.gradeMmS] balance quality grade G [mm/s] (6.3 for G 6.3)
 * @param {number} [o.eccentricity] mass eccentricity e [m] (alternative to the grade)
 * @param {number} [o.omega] rotational speed [rad/s], or
 * @param {number} [o.speedRpm] rotational speed [rev/min]
 * @param {number} [o.serviceFactor=1] multiplier on the force [-]
 * @returns {{force:number, eccentricity:number, unbalance:number, omega:number, frequencyHz:number}} force amplitude [N], eccentricity [m],
 *   unbalance m·e [kg·m], speed [rad/s], rotation frequency [Hz]
 */
export function unbalanceForce({ mass, gradeMmS, eccentricity, omega, speedRpm, serviceFactor = 1 }) {
  pos('mass', mass); pos('serviceFactor', serviceFactor);
  const w = speedOf({ omega, speedRpm });
  let e;
  if (eccentricity !== undefined) e = nonNeg('eccentricity', eccentricity);
  else e = pos('gradeMmS', gradeMmS) * 1e-3 / w;
  return { force: serviceFactor * mass * e * w * w, eccentricity: e, unbalance: mass * e, omega: w, frequencyHz: w / (2 * Math.PI) };
}

/* ------------------------------------------------------------------ SDOF response and isolation */

/**
 * Dynamic magnification and force transmissibility of a viscously damped SDOF system at frequency ratio r = ω/ωn:
 * D = 1/√[(1 − r²)² + (2ζr)²], TR = √[1 + (2ζr)²]·D (Den Hartog ch. 2). Isolation (TR < 1) requires r > √2.
 * @param {object} o
 * @param {number} o.r frequency ratio ω/ωn [-]
 * @param {number} [o.zeta=0] damping ratio [-]
 * @returns {{magnification:number, transmissibility:number, isolationEfficiency:number, phase:number}} D [-], TR [-], 1 − TR [-] (negative in
 *   the amplification range), phase lag of the displacement behind the force [rad]
 */
export function transmissibility({ r, zeta = 0 }) {
  nonNeg('r', r); nonNeg('zeta', zeta);
  const den = Math.sqrt((1 - r * r) ** 2 + (2 * zeta * r) ** 2);
  const D = den === 0 ? Infinity : 1 / den;
  const TR = Math.sqrt(1 + (2 * zeta * r) ** 2) * D;
  return { magnification: D, transmissibility: TR, isolationEfficiency: 1 - TR, phase: Math.atan2(2 * zeta * r, 1 - r * r) };
}

/**
 * Steady-state response of a machine–foundation SDOF system to a harmonic force: m ẍ + c ẋ + k x = F0 sin ωt.
 * For a rotating unbalance F0 = m_u e ω² (pass `unbalance` = m_u·e) the amplitude is X = (m_u e/m) r² D.
 * @param {object} o
 * @param {number} o.mass total vibrating mass (machine + block) [kg]
 * @param {number} o.k support stiffness [N/m]
 * @param {number} [o.zeta] damping ratio [-], or
 * @param {number} [o.c] viscous damping coefficient [N·s/m]
 * @param {number} [o.omega] forcing frequency [rad/s] (or speedRpm / frequencyHz)
 * @param {number} [o.speedRpm] forcing frequency as a speed [rev/min]
 * @param {number} [o.frequencyHz] forcing frequency [Hz]
 * @param {number} [o.F0] force amplitude [N], or
 * @param {number} [o.unbalance] rotating unbalance m_u·e [kg·m]
 * @returns {{amplitude:number, velocity:number, velocityRms:number, acceleration:number, phase:number, r:number, magnification:number,
 *   transmissibility:number, transmittedForce:number, isolationEfficiency:number, omegaN:number, naturalHz:number, zeta:number, F0:number,
 *   staticDeflection:number}} displacement amplitude [m], velocity amplitude and rms [m/s], acceleration amplitude [m/s²], phase lag [rad],
 *   frequency ratio, D, TR, force transmitted to the support [N], 1 − TR, natural frequency [rad/s, Hz], damping ratio, force amplitude [N],
 *   static deflection under the weight m g/k [m]
 */
export function sdofForcedResponse({ mass, k, zeta, c, omega, speedRpm, frequencyHz, F0, unbalance }) {
  pos('mass', mass); pos('k', k);
  const w = speedOf({ omega, speedRpm, frequencyHz });
  const wn = Math.sqrt(k / mass);
  let z = 0;
  if (zeta !== undefined) z = nonNeg('zeta', zeta); else if (c !== undefined) z = nonNeg('c', c) / (2 * mass * wn);
  let F;
  if (F0 !== undefined) F = nonNeg('F0', F0); else if (unbalance !== undefined) F = nonNeg('unbalance', unbalance) * w * w; else throw new RangeError('F0 or unbalance is required');
  const t = transmissibility({ r: w / wn, zeta: z });
  const X = F / k * t.magnification;
  return {
    amplitude: X, velocity: X * w, velocityRms: X * w / Math.SQRT2, acceleration: X * w * w, phase: t.phase, r: w / wn, magnification: t.magnification,
    transmissibility: t.transmissibility, transmittedForce: F * t.transmissibility, isolationEfficiency: t.isolationEfficiency, omegaN: wn,
    naturalHz: wn / (2 * Math.PI), zeta: z, F0: F, staticDeflection: mass * G0 / k,
  };
}

/**
 * Isolator selection for a required isolation efficiency (undamped mounts): TR = 1/(r² − 1) for r > √2, so
 * r = √(1 + 1/TR), fn = f/r, and the static deflection of linear mounts under the machine weight δst = g/(2π fn)².
 * Damping raises the transmissibility in the isolation range; check the result with transmissibility().
 * @param {object} o
 * @param {number} o.efficiency required isolation efficiency 1 − TR [-], in (0, 1)
 * @param {number} [o.frequencyHz] lowest disturbing frequency [Hz] (or speedRpm / omega)
 * @param {number} [o.speedRpm] [rev/min]
 * @param {number} [o.omega] [rad/s]
 * @param {number} [o.mass] supported mass [kg] for the required total stiffness
 * @returns {{naturalHz:number, r:number, staticDeflection:number, stiffness:(number|undefined)}} required natural frequency of the mounted
 *   system [Hz] (upper limit), frequency ratio [-], static deflection [m], total mount stiffness [N/m] when the mass is given
 */
export function isolatorSelection({ efficiency, frequencyHz, speedRpm, omega, mass }) {
  if (!(efficiency > 0 && efficiency < 1)) throw new RangeError('efficiency must lie in (0, 1)');
  const w = speedOf({ omega, speedRpm, frequencyHz });
  const r = Math.sqrt(1 + 1 / (1 - efficiency));
  const wn = w / r;
  return { naturalHz: wn / (2 * Math.PI), r, staticDeflection: G0 / (wn * wn), stiffness: mass === undefined ? undefined : pos('mass', mass) * wn * wn };
}

/* ------------------------------------------------------------------ 2-DOF machine / foundation */

/**
 * Two-degree-of-freedom model of a machine (m1) on isolators (k1, c1) carried by a foundation block or deck
 * structure (m2) on its own support (k2, c2), with a harmonic force F0 on the machine:
 *   m1 ẍ1 + c1(ẋ1 − ẋ2) + k1(x1 − x2) = F0 e^(iωt),   m2 ẍ2 + c2 ẋ2 + k2 x2 − c1(ẋ1 − ẋ2) − k1(x1 − x2) = 0.
 * Solved in complex form; the undamped natural frequencies are the roots of
 * m1 m2 λ² − [m1 (k1 + k2) + m2 k1] λ + k1 k2 = 0, λ = ω² (Den Hartog ch. 3).
 * @param {object} o
 * @param {number} o.m1 machine mass [kg]
 * @param {number} o.k1 isolator stiffness [N/m]
 * @param {number} [o.c1=0] isolator damping [N·s/m]
 * @param {number} o.m2 foundation / supporting structure mass [kg]
 * @param {number} o.k2 foundation support stiffness [N/m]
 * @param {number} [o.c2=0] foundation damping [N·s/m]
 * @param {number} o.F0 force amplitude on the machine [N]
 * @param {number} [o.omega] forcing frequency [rad/s] (or speedRpm / frequencyHz)
 * @param {number} [o.speedRpm] [rev/min]
 * @param {number} [o.frequencyHz] [Hz]
 * @returns {{x1:number, x2:number, phase1:number, phase2:number, relative:number, v1:number, v2:number, forceToGround:number,
 *   forceInIsolators:number, transmissibility:number, naturalOmegas:[number, number], naturalHz:[number, number]}} displacement amplitudes
 *   of machine and foundation [m], phase lags [rad], amplitude of the isolator deformation [m], velocity amplitudes [m/s], force amplitudes
 *   transmitted to the ground and through the isolators [N], forceToGround/F0 [-], undamped natural frequencies [rad/s] and [Hz]
 */
export function twoDofFoundation({ m1, k1, c1 = 0, m2, k2, c2 = 0, F0, omega, speedRpm, frequencyHz }) {
  pos('m1', m1); pos('k1', k1); nonNeg('c1', c1); pos('m2', m2); pos('k2', k2); nonNeg('c2', c2); nonNeg('F0', F0);
  const w = speedOf({ omega, speedRpm, frequencyHz });
  // complex 2×2: [a, −b; −b, d] {X1, X2} = {F0, 0},  a = k1 − m1ω² + iωc1, b = k1 + iωc1, d = k1 + k2 − m2ω² + iω(c1 + c2)
  const ar = k1 - m1 * w * w; const ai = w * c1;
  const br = k1; const bi = w * c1;
  const dr = k1 + k2 - m2 * w * w; const di = w * (c1 + c2);
  const detR = ar * dr - ai * di - (br * br - bi * bi);
  const detI = ar * di + ai * dr - 2 * br * bi;
  const den = detR * detR + detI * detI;
  if (den === 0) throw new RangeError('undamped system excited exactly at a natural frequency: response is unbounded');
  const div = (xr, xi) => [(xr * detR + xi * detI) / den, (xi * detR - xr * detI) / den];
  const [x1r, x1i] = div(F0 * dr, F0 * di);
  const [x2r, x2i] = div(F0 * br, F0 * bi);
  const x1 = Math.hypot(x1r, x1i); const x2 = Math.hypot(x2r, x2i);
  const rel = Math.hypot(x1r - x2r, x1i - x2i);
  const B = m1 * (k1 + k2) + m2 * k1;
  const disc = Math.sqrt(B * B - 4 * m1 * m2 * k1 * k2);
  const w1 = Math.sqrt((B - disc) / (2 * m1 * m2)); const w2 = Math.sqrt((B + disc) / (2 * m1 * m2));
  const fg = Math.hypot(k2, w * c2) * x2;
  return {
    x1, x2, phase1: Math.atan2(-x1i, x1r), phase2: Math.atan2(-x2i, x2r), relative: rel, v1: x1 * w, v2: x2 * w, forceToGround: fg,
    forceInIsolators: Math.hypot(k1, w * c1) * rel, transmissibility: F0 > 0 ? fg / F0 : 0, naturalOmegas: [w1, w2],
    naturalHz: [w1 / (2 * Math.PI), w2 / (2 * Math.PI)],
  };
}

/* ------------------------------------------------------------------ footing on a half-space */

/**
 * Equivalent radii of a rectangular footing B × L for the circular-footing formulas (Richart, Hall & Woods 1970, ch. 10):
 * translation r0 = √(B L/π); rocking about the axis parallel to B, r0 = [B L³/(3π)]^¼; torsion r0 = [B L (B² + L²)/(6π)]^¼.
 * @param {object} o
 * @param {number} o.B footing dimension parallel to the rocking axis [m]
 * @param {number} o.L footing dimension perpendicular to the rocking axis [m]
 * @returns {{translation:number, rocking:number, torsion:number}} radii [m]
 */
export function equivalentFootingRadius({ B, L }) {
  pos('B', B); pos('L', L);
  return { translation: Math.sqrt(B * L / Math.PI), rocking: (B * L ** 3 / (3 * Math.PI)) ** 0.25, torsion: (B * L * (B * B + L * L) / (6 * Math.PI)) ** 0.25 };
}

/**
 * Frequency-independent spring constants and equivalent damping ratios of a rigid circular footing on an elastic
 * half-space — the lumped-parameter analogue of Richart, Hall & Woods (1970), ch. 10:
 *   vertical    kz = 4 G r/(1 − ν),         Bz = (1 − ν) m/(4 ρ r³),             Dz = 0.425/√Bz
 *   horizontal  kx = 8 G r/(2 − ν),         Bx = (7 − 8ν) m/(32 (1 − ν) ρ r³),   Dx = 0.288/√Bx
 *   rocking     kψ = 8 G r³/(3 (1 − ν)),    Bψ = 3 (1 − ν) Iψ/(8 ρ r⁵),          Dψ = 0.15/[(1 + Bψ) √Bψ]
 *   torsion     kθ = 16 G r³/3,             Bθ = Iθ/(ρ r⁵),                      Dθ = 0.5/(1 + 2 Bθ)
 * (B: mass / inertia ratio, D: geometric damping ratio). Surface footing; embedment and soil material damping are not
 * included. Use equivalentFootingRadius for rectangular bases (a separate radius for each mode).
 * @param {object} o
 * @param {number} o.G soil shear modulus [Pa]
 * @param {number} o.nu soil Poisson's ratio [-]
 * @param {number} o.rho soil mass density [kg/m³]
 * @param {number} o.radius footing radius [m] (used for all modes unless the mode-specific radii are given)
 * @param {number} [o.radiusRocking] equivalent radius for rocking [m]
 * @param {number} [o.radiusTorsion] equivalent radius for torsion [m]
 * @param {number} [o.mass] machine + foundation mass [kg] (for the translational damping ratios)
 * @param {number} [o.inertiaRocking] mass moment of inertia about the rocking axis through the base [kg·m²]
 * @param {number} [o.inertiaTorsion] mass moment of inertia about the vertical axis [kg·m²]
 * @returns {{kz:number, kx:number, kRocking:number, kTorsion:number, Bz:(number|undefined), Bx:(number|undefined), BRocking:(number|undefined),
 *   BTorsion:(number|undefined), Dz:(number|undefined), Dx:(number|undefined), DRocking:(number|undefined), DTorsion:(number|undefined)}}
 *   stiffnesses [N/m, N·m/rad], mass ratios [-], damping ratios [-]
 */
export function footingImpedance({ G, nu, rho, radius, radiusRocking, radiusTorsion, mass, inertiaRocking, inertiaTorsion }) {
  pos('G', G); pos('rho', rho); pos('radius', radius);
  if (!(nu >= 0 && nu < 0.5 + 1e-12)) throw new RangeError('nu must lie in [0, 0.5]');
  const rr = radiusRocking === undefined ? radius : pos('radiusRocking', radiusRocking);
  const rt = radiusTorsion === undefined ? radius : pos('radiusTorsion', radiusTorsion);
  const out = {
    kz: 4 * G * radius / (1 - nu), kx: 8 * G * radius / (2 - nu), kRocking: 8 * G * rr ** 3 / (3 * (1 - nu)), kTorsion: 16 * G * rt ** 3 / 3,
    Bz: undefined, Bx: undefined, BRocking: undefined, BTorsion: undefined, Dz: undefined, Dx: undefined, DRocking: undefined, DTorsion: undefined,
  };
  if (mass !== undefined) {
    pos('mass', mass);
    out.Bz = (1 - nu) / 4 * mass / (rho * radius ** 3); out.Dz = 0.425 / Math.sqrt(out.Bz);
    out.Bx = (7 - 8 * nu) / (32 * (1 - nu)) * mass / (rho * radius ** 3); out.Dx = 0.288 / Math.sqrt(out.Bx);
  }
  if (inertiaRocking !== undefined) {
    pos('inertiaRocking', inertiaRocking);
    out.BRocking = 3 * (1 - nu) / 8 * inertiaRocking / (rho * rr ** 5); out.DRocking = 0.15 / ((1 + out.BRocking) * Math.sqrt(out.BRocking));
  }
  if (inertiaTorsion !== undefined) {
    pos('inertiaTorsion', inertiaTorsion);
    out.BTorsion = inertiaTorsion / (rho * rt ** 5); out.DTorsion = 0.5 / (1 + 2 * out.BTorsion);
  }
  return out;
}

/* ------------------------------------------------------------------ acceptance checks */

/**
 * Separation margin between an excitation frequency and a natural frequency: margin = |fe − fn|/fn.
 * The required margin is a project / code criterion supplied by the caller (±20 % is a widely used value for machine
 * foundations and supporting structures); `tuning` reports whether the support is low-tuned (fn < fe) or high-tuned.
 * @param {object} o
 * @param {number} o.excitationHz excitation frequency [Hz] (running speed or a multiple)
 * @param {number} o.naturalHz natural frequency [Hz]
 * @param {number} [o.required=0.2] required margin [-]
 * @returns {{margin:number, ratio:number, ok:boolean, tuning:'low'|'high', avoidBandHz:[number, number]}} margin [-], fe/fn [-], pass flag,
 *   tuning, and the band of natural frequencies to avoid for this excitation [Hz]
 */
export function separationMargin({ excitationHz, naturalHz, required = 0.2 }) {
  pos('excitationHz', excitationHz); pos('naturalHz', naturalHz); nonNeg('required', required);
  const margin = Math.abs(excitationHz - naturalHz) / naturalHz;
  // |fe − fn|/fn < required  ⇔  fe/(1 + required) < fn < fe/(1 − required)
  return {
    margin, ratio: excitationHz / naturalHz, ok: margin >= required, tuning: naturalHz < excitationHz ? 'low' : 'high',
    avoidBandHz: [excitationHz / (1 + required), required < 1 ? excitationHz / (1 - required) : Infinity],
  };
}

/**
 * Vibration severity zone from the broadband rms velocity and the zone boundaries of the applicable standard
 * (ISO 20816 / ISO 10816 series: A — newly commissioned, B — acceptable for unrestricted long-term operation,
 * C — unsatisfactory for long-term operation, D — may cause damage). The boundaries depend on machine type, size and
 * support flexibility and must be supplied from the relevant part of the standard or the project specification.
 * @param {object} o
 * @param {number} o.velocityRms measured or predicted rms vibration velocity [m/s]
 * @param {{AB:number, BC:number, CD:number}} o.limits zone boundaries [m/s], ascending
 * @returns {{zone:'A'|'B'|'C'|'D', utilisation:number, acceptable:boolean}} zone, velocity relative to the B/C boundary [-], zone A or B
 */
export function vibrationSeverity({ velocityRms, limits }) {
  nonNeg('velocityRms', velocityRms);
  if (!limits || !(limits.AB > 0 && limits.BC > limits.AB && limits.CD > limits.BC)) throw new RangeError('limits must be { AB, BC, CD } in ascending order');
  const zone = velocityRms <= limits.AB ? 'A' : (velocityRms <= limits.BC ? 'B' : (velocityRms <= limits.CD ? 'C' : 'D'));
  return { zone, utilisation: velocityRms / limits.BC, acceptable: zone === 'A' || zone === 'B' };
}

/**
 * Conversion between the amplitudes of harmonic displacement, velocity and acceleration at one frequency:
 * v = ω x, a = ω² x; rms = peak/√2; peak-to-peak displacement = 2 x.
 * @param {object} o
 * @param {number} o.frequencyHz frequency [Hz]
 * @param {number} [o.displacement] displacement amplitude (zero-to-peak) [m], or
 * @param {number} [o.velocity] velocity amplitude [m/s], or
 * @param {number} [o.acceleration] acceleration amplitude [m/s²]
 * @returns {{displacement:number, displacementPeakToPeak:number, velocity:number, velocityRms:number, acceleration:number, accelerationRms:number}}
 */
export function harmonicAmplitudes({ frequencyHz, displacement, velocity, acceleration }) {
  const w = 2 * Math.PI * pos('frequencyHz', frequencyHz);
  let x;
  if (displacement !== undefined) x = nonNeg('displacement', displacement);
  else if (velocity !== undefined) x = nonNeg('velocity', velocity) / w;
  else if (acceleration !== undefined) x = nonNeg('acceleration', acceleration) / (w * w);
  else throw new RangeError('one of displacement, velocity or acceleration is required');
  return {
    displacement: x, displacementPeakToPeak: 2 * x, velocity: x * w, velocityRms: x * w / Math.SQRT2, acceleration: x * w * w,
    accelerationRms: x * w * w / Math.SQRT2,
  };
}

/* ------------------------------------------------------------------ Campbell diagram */

/**
 * Data for a Campbell (interference) diagram: natural frequencies against rotational speed, excitation order lines
 * f = n·N/60 and their intersections (potential resonances). Natural frequencies may be constants or functions of the
 * speed (gyroscopic splitting of whirl modes, speed-dependent bearing or centrifugal stiffening); crossings are located
 * by root finding between the sampled speeds.
 * @param {object} o
 * @param {Array<number|function(number):number>} o.modes natural frequencies [Hz], constants or functions of the speed [rev/min]
 * @param {number[]} [o.orders=[1, 2]] excitation orders (multiples of the running speed; e.g. blade-pass = number of blades)
 * @param {number} o.speedMinRpm lowest speed of the diagram [rev/min] (>= 0)
 * @param {number} o.speedMaxRpm highest speed [rev/min]
 * @param {[number, number]} [o.operatingRangeRpm] operating speed range for the margin check [rev/min]
 * @param {number} [o.required=0.2] required separation margin within the operating range [-] (project criterion)
 * @param {number} [o.nPoints=101] number of speeds sampled
 * @returns {{speedRpm:Float64Array, modeHz:Float64Array[], orderHz:Float64Array[], crossings:Array<{mode:number, order:number, speedRpm:number,
 *   hz:number, inOperatingRange:boolean}>, margins:Array<{mode:number, order:number, margin:number, speedRpm:number, ok:boolean}>}}
 *   sampled speeds, natural frequency of each mode and frequency of each order line at those speeds [Hz], crossings (critical speeds),
 *   and — when an operating range is given — the smallest separation margin |f_exc − f_n|/f_n of each mode/order pair within the range
 */
export function campbellData({ modes, orders = [1, 2], speedMinRpm, speedMaxRpm, operatingRangeRpm, required = 0.2, nPoints = 101 }) {
  if (!Array.isArray(modes) || modes.length === 0) throw new RangeError('modes must be a non-empty array');
  if (!Array.isArray(orders) || orders.length === 0) throw new RangeError('orders must be a non-empty array');
  nonNeg('speedMinRpm', speedMinRpm); pos('speedMaxRpm', speedMaxRpm);
  if (!(speedMaxRpm > speedMinRpm)) throw new RangeError('speedMaxRpm must exceed speedMinRpm');
  if (!Number.isInteger(nPoints) || nPoints < 2) throw new RangeError('nPoints must be an integer >= 2');
  const fn = modes.map((m, i) => {
    if (typeof m === 'function') return m;
    pos(`modes[${i}]`, m);
    return () => m;
  });
  orders.forEach((n, i) => pos(`orders[${i}]`, n));
  const speed = new Float64Array(nPoints);
  for (let i = 0; i < nPoints; i++) speed[i] = speedMinRpm + (speedMaxRpm - speedMinRpm) * i / (nPoints - 1);
  const modeHz = fn.map((f) => Float64Array.from(speed, (s) => f(s)));
  const orderHz = orders.map((n) => Float64Array.from(speed, (s) => n * s / 60));
  const inRange = (s) => !!operatingRangeRpm && s >= operatingRangeRpm[0] && s <= operatingRangeRpm[1];
  const crossings = [];
  fn.forEach((f, mi) => {
    orders.forEach((n) => {
      const g = (s) => n * s / 60 - f(s);
      for (let i = 0; i + 1 < nPoints; i++) {
        const a = g(speed[i]); const b = g(speed[i + 1]);
        if (a === 0) crossings.push({ mode: mi, order: n, speedRpm: speed[i], hz: f(speed[i]), inOperatingRange: inRange(speed[i]) });
        else if (a * b < 0) {
          const s = brent(g, speed[i], speed[i + 1], { tol: 1e-12 });
          crossings.push({ mode: mi, order: n, speedRpm: s, hz: f(s), inOperatingRange: inRange(s) });
        }
      }
      if (g(speed[nPoints - 1]) === 0) crossings.push({ mode: mi, order: n, speedRpm: speed[nPoints - 1], hz: f(speed[nPoints - 1]), inOperatingRange: inRange(speed[nPoints - 1]) });
    });
  });
  const margins = [];
  if (operatingRangeRpm) {
    const [lo, hi] = operatingRangeRpm;
    if (!(hi >= lo && lo >= 0)) throw new RangeError('operatingRangeRpm must be [min, max]');
    fn.forEach((f, mi) => {
      orders.forEach((n) => {
        let best = Infinity; let at = lo;
        const m = 200;
        for (let i = 0; i <= m; i++) {
          const s = lo + (hi - lo) * i / m;
          const fnat = f(s);
          const mg = Math.abs(n * s / 60 - fnat) / fnat;
          if (mg < best) { best = mg; at = s; }
        }
        if (crossings.some((c) => c.mode === mi && c.order === n && c.inOperatingRange)) best = 0;
        margins.push({ mode: mi, order: n, margin: best, speedRpm: at, ok: best >= required });
      });
    });
  }
  return { speedRpm: speed, modeHz, orderHz, crossings, margins };
}
