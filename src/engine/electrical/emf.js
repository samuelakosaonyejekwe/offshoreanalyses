/**
 * Electromagnetic fields and induction around submarine cables.
 *
 * References
 *  - Biot–Savart law for long straight conductors, B = μ0·I/(2π·r) (any electromagnetics text).
 *  - ICNIRP (2010) Guidelines for limiting exposure to time-varying electric and magnetic fields (1 Hz–100 kHz):
 *    general-public reference level 200 µT at 50 Hz (100 µT in the 1998 guidelines); exposed as arguments.
 *  - Carson (1926), Bell Syst. Tech. J. 5: mutual impedance of parallel conductors with earth return,
 *    first-term approximation Zm = ωμ0/8 + jωμ0/(2π)·ln(De/d), De = 658.87·sqrt(ρ/f).
 *  - EN 15280 / ISO 18086 (AC corrosion likelihood of buried pipelines; current-density criteria as arguments).
 *
 * Geometry: 2-D cross-section, x horizontal (m), z vertical (m, positive up, seabed at z = 0, so a buried cable
 * has negative z). Conductors run along y. Phasor currents are rms complex { re, im }; a current is positive
 * along +y. SI units, flux density in tesla.
 */

import { cx, cadd, cmul, cscale, cabs, cpolar, toCx } from './complex.js';

const MU0 = 4e-7 * Math.PI;
const TWO_PI = 2 * Math.PI;

function pos(v, name) {
  if (!(v > 0) || !Number.isFinite(v)) throw new RangeError(`${name}: must be a positive finite number`);
}
function nonNeg(v, name) {
  if (!(v >= 0) || !Number.isFinite(v)) throw new RangeError(`${name}: must be a non-negative finite number`);
}

/**
 * Magnetic flux density at a point from a set of long, straight, parallel conductors (Biot–Savart superposition).
 * For a current I along +y at (xc, zc): Bx = μ0·I/(2π)·(z − zc)/r², Bz = −μ0·I/(2π)·(x − xc)/r².
 * With phasor currents the field is elliptically polarised; `rms` = sqrt(|Bx|² + |Bz|²) is the rms resultant
 * (the quantity compared with exposure reference levels) and `major` the rms value along the major axis.
 * For DC give real currents: `rms` is then the static field magnitude.
 * @param {object} o
 * @param {{x:number, z:number, I:number|{re:number, im:number}}[]} o.conductors position (m) and current (A)
 * @param {{x:number, z:number}} o.point field point, m
 * @param {number} [o.shieldingFactor=1] ratio of shielded to unshielded field (≤ 1) for sheath/armour screening
 *        (steel-wire armour and lead sheath typically give 0.8–1.0 at 50 Hz for the far field of a three-core
 *        cable; measured or FEA value)
 * @returns {{Bx:{re:number, im:number}, Bz:{re:number, im:number}, rms:number, major:number, minor:number}} tesla
 */
export function fluxDensity({ conductors, point, shieldingFactor = 1 }) {
  if (!Array.isArray(conductors) || conductors.length === 0) throw new RangeError('conductors: must be a non-empty array');
  if (!(shieldingFactor >= 0 && shieldingFactor <= 1)) throw new RangeError('shieldingFactor: must be in [0, 1]');
  let Bx = cx(0, 0); let Bz = cx(0, 0);
  for (const c of conductors) {
    const dx = point.x - c.x; const dz = point.z - c.z;
    const r2 = dx * dx + dz * dz;
    if (!(r2 > 0)) throw new RangeError('point: coincides with a conductor axis');
    const k = (MU0 / TWO_PI) / r2;
    const I = toCx(c.I);
    Bx = cadd(Bx, cscale(I, k * dz));
    Bz = cadd(Bz, cscale(I, -k * dx));
  }
  Bx = cscale(Bx, shieldingFactor); Bz = cscale(Bz, shieldingFactor);
  const sumSq = Bx.re * Bx.re + Bx.im * Bx.im + Bz.re * Bz.re + Bz.im * Bz.im;
  // Semi-axes of the field ellipse from the complex "square" Bx² + Bz².
  const sq = cadd(cmul(Bx, Bx), cmul(Bz, Bz));
  const m = cabs(sq);
  return { Bx, Bz, rms: Math.sqrt(sumSq), major: Math.sqrt(0.5 * (sumSq + m)), minor: Math.sqrt(Math.max(0, 0.5 * (sumSq - m))) };
}

/**
 * Conductor sets for standard arrangements, centred on x = 0 at height z (negative = buried).
 *  - 'flat'     three phases at x = −s, 0, +s, currents I∠0, I∠−120°, I∠+120°
 *  - 'trefoil'  three phases at the corners of an equilateral triangle of side s (apex up), centroid at z
 *  - 'bipole'   two DC conductors at x = ∓s/2 with currents +I and −I (bundled bipole: s = cable diameter)
 *  - 'monopole' a single DC conductor (sea or remote return)
 * @param {object} o
 * @param {'flat'|'trefoil'|'bipole'|'monopole'} o.arrangement
 * @param {number} o.I rms phase current (AC) or pole current (DC), A
 * @param {number} [o.s=0] axial spacing, m
 * @param {number} [o.z=0] height of the group centre relative to the seabed, m
 * @returns {{x:number, z:number, I:{re:number, im:number}}[]}
 */
export function arrangement({ arrangement: kind, I, s = 0, z = 0 }) {
  if (!Number.isFinite(I)) throw new RangeError('I: must be finite');
  nonNeg(s, 's');
  const a = (ang) => cpolar(I, ang);
  const d120 = TWO_PI / 3;
  switch (kind) {
    case 'flat':
      pos(s, 's');
      return [{ x: -s, z, I: a(0) }, { x: 0, z, I: a(-d120) }, { x: s, z, I: a(d120) }];
    case 'trefoil': {
      pos(s, 's');
      const h = s / Math.sqrt(3);
      return [{ x: 0, z: z + h, I: a(0) }, { x: -s / 2, z: z - h / 2, I: a(-d120) }, { x: s / 2, z: z - h / 2, I: a(d120) }];
    }
    case 'bipole':
      pos(s, 's');
      return [{ x: -s / 2, z, I: cx(I, 0) }, { x: s / 2, z, I: cx(-I, 0) }];
    case 'monopole':
      return [{ x: 0, z, I: cx(I, 0) }];
    default:
      throw new RangeError('arrangement: must be flat, trefoil, bipole or monopole');
  }
}

/**
 * Indicative attenuation of the external field of a three-core cable by the helical lay of its cores, as a
 * multiplying factor on the field of the equivalent untwisted trefoil.
 * The field of a twisted triplet is governed by modified Bessel functions K(2π·r/p) of the radial distance r and
 * lay length p (Pettersson & Schönborg, IEEE Trans. Power Delivery 12 (1997)), whose asymptotic behaviour is
 * exp(−2π·r/p). Only that exponential term is evaluated here, as min(1, exp(−2π·(r − a)/p)) with a the radius of
 * the core centres; the algebraic prefactors are omitted. Treat the result as a screening-level factor and
 * replace it by a measured or computed value (shieldingFactor argument of fluxDensity) for assessments.
 * @param {object} o
 * @param {number} o.r distance from the cable axis, m
 * @param {number} o.pitch lay length of the cores, m
 * @param {number} [o.coreRadius=0] radius of the core centres from the cable axis a, m
 * @returns {number} factor (0–1)
 */
export function twistAttenuation({ r, pitch, coreRadius = 0 }) {
  pos(r, 'r'); pos(pitch, 'pitch'); nonNeg(coreRadius, 'coreRadius');
  return Math.min(1, Math.exp((-TWO_PI * Math.max(0, r - coreRadius)) / pitch));
}

/**
 * Flux-density profile along a horizontal line at a height above the seabed.
 * @param {object} o
 * @param {{x:number, z:number, I:number|{re:number, im:number}}[]} o.conductors (buried conductors have z < 0)
 * @param {number} [o.height=0] height of the profile above the seabed, m
 * @param {number} [o.xMin=-20] @param {number} [o.xMax=20] m @param {number} [o.n=201] points
 * @param {number} [o.shieldingFactor=1]
 * @returns {{x:Float64Array, B:Float64Array, Bmax:number, xMax:number}} rms resultant flux density (T)
 */
export function fluxDensityProfile({ conductors, height = 0, xMin = -20, xMax = 20, n = 201, shieldingFactor = 1 }) {
  if (!(n >= 2)) throw new RangeError('n: at least 2 points');
  if (!(xMax > xMin)) throw new RangeError('xMax: must exceed xMin');
  const x = new Float64Array(n); const B = new Float64Array(n);
  let Bmax = 0; let xAt = 0;
  for (let i = 0; i < n; i++) {
    x[i] = xMin + ((xMax - xMin) * i) / (n - 1);
    B[i] = fluxDensity({ conductors, point: { x: x[i], z: height }, shieldingFactor }).rms;
    if (B[i] > Bmax) { Bmax = B[i]; xAt = x[i]; }
  }
  return { x, B, Bmax, xMax: xAt };
}

/**
 * Distance from the cable group at which the flux density falls below a threshold, searched along a ray from
 * the group centre (default vertically upwards). Bisection on the monotonic far part of the field.
 * Typical thresholds: 200 µT (ICNIRP 2010, general public, 50 Hz), 100 µT (ICNIRP 1998), or a project value for
 * marine fauna; geomagnetic background is about 25–65 µT.
 * @param {object} o
 * @param {{x:number, z:number, I:number|{re:number, im:number}}[]} o.conductors
 * @param {number} o.threshold flux density, T
 * @param {number} [o.directionDeg=90] direction of the ray from the +x axis towards +z, degrees
 * @param {number} [o.shieldingFactor=1] @param {number} [o.rMax=1000] search limit, m
 * @returns {number} distance from the centroid of the conductors, m (0 if the field is below the threshold
 *          everywhere outside the group envelope)
 */
export function distanceToThreshold({ conductors, threshold, directionDeg = 90, shieldingFactor = 1, rMax = 1000 }) {
  pos(threshold, 'threshold');
  let cxm = 0; let czm = 0; let env = 0;
  for (const c of conductors) { cxm += c.x; czm += c.z; }
  cxm /= conductors.length; czm /= conductors.length;
  for (const c of conductors) env = Math.max(env, Math.hypot(c.x - cxm, c.z - czm));
  const ux = Math.cos((directionDeg * Math.PI) / 180); const uz = Math.sin((directionDeg * Math.PI) / 180);
  const B = (r) => fluxDensity({ conductors, point: { x: cxm + r * ux, z: czm + r * uz }, shieldingFactor }).rms;
  let lo = env * 1.0001 + 1e-6;
  if (B(lo) <= threshold) return 0;
  if (B(rMax) > threshold) throw new RangeError('rMax: field still above the threshold at the search limit');
  let hi = rMax;
  for (let i = 0; i < 200 && hi - lo > 1e-9 * hi; i++) {
    const mid = Math.sqrt(lo * hi);
    if (B(mid) > threshold) lo = mid; else hi = mid;
  }
  return 0.5 * (lo + hi);
}

/**
 * Static field of DC cables superposed on the geomagnetic field, and the resulting compass deviation.
 * The cable runs horizontally with heading ψ (clockwise from magnetic north); the cross-section x axis points to
 * the right of the heading and z is up. Geomagnetic field: total intensity F, inclination I (positive down);
 * horizontal component H = F·cos I towards magnetic north, vertical component F·sin I downwards.
 * Compass deviation is the rotation of the horizontal field direction.
 * @param {object} o
 * @param {{x:number, z:number, I:number}[]} o.conductors DC conductors (real currents, A)
 * @param {{x:number, z:number}} o.point field point, m
 * @param {number} o.headingDeg cable heading from magnetic north, degrees
 * @param {number} [o.geomagneticTotal=50e-6] F, T (about 25–65 µT depending on location)
 * @param {number} [o.inclinationDeg=65] magnetic inclination (dip), degrees, positive downwards
 * @returns {{cable:number, total:number, deltaTotal:number, compassDeviationDeg:number, horizontal:number}}
 *          cable field magnitude, total field magnitude, change of total intensity (T), deviation of the
 *          horizontal field direction (degrees, positive clockwise/east) and resulting horizontal intensity (T)
 */
export function dcFieldWithGeomagnetic({ conductors, point, headingDeg, geomagneticTotal = 50e-6, inclinationDeg = 65 }) {
  nonNeg(geomagneticTotal, 'geomagneticTotal');
  if (!Number.isFinite(headingDeg)) throw new RangeError('headingDeg: must be finite');
  const b = fluxDensity({ conductors, point });
  const bx = b.Bx.re; const bz = b.Bz.re;
  const psi = (headingDeg * Math.PI) / 180; const inc = (inclinationDeg * Math.PI) / 180;
  // The x axis points to the right of the heading, i.e. towards azimuth ψ + 90°.
  const bn = bx * Math.cos(psi + Math.PI / 2); const be = bx * Math.sin(psi + Math.PI / 2);
  const H = geomagneticTotal * Math.cos(inc); const Zdown = geomagneticTotal * Math.sin(inc);
  const n = H + bn; const e = be; const up = -Zdown + bz;
  const total = Math.hypot(n, e, up);
  return { cable: Math.hypot(bx, bz), total, deltaTotal: total - geomagneticTotal,
    compassDeviationDeg: (Math.atan2(e, n) * 180) / Math.PI, horizontal: Math.hypot(n, e) };
}

/**
 * Motionally induced electric field in seawater or in an animal moving through a magnetic field, E = v × B;
 * magnitude E = v·B·sin(angle).
 * @param {object} o
 * @param {number} o.velocity speed of the water or animal relative to the field, m/s
 * @param {number} o.B flux density, T
 * @param {number} [o.angleDeg=90] angle between velocity and field, degrees
 * @returns {number} electric field, V/m
 */
export function motionalElectricField({ velocity, B, angleDeg = 90 }) {
  nonNeg(velocity, 'velocity'); nonNeg(B, 'B');
  return velocity * B * Math.abs(Math.sin((angleDeg * Math.PI) / 180));
}

/**
 * Electric field induced by a uniform alternating flux density in a conducting loop or body of circular
 * cross-section (Faraday's law, E·2πa = πa²·dB/dt): E = π·f·a·B, with E and B both rms (or both peak).
 * @param {object} o
 * @param {number} o.B flux density, T @param {number} o.f frequency, Hz
 * @param {number} o.radius radius of the induction loop a, m (e.g. half the body width of an animal)
 * @returns {number} induced electric field at the loop radius, V/m
 */
export function faradayElectricField({ B, f, radius }) {
  nonNeg(B, 'B'); nonNeg(f, 'f'); nonNeg(radius, 'radius');
  return Math.PI * f * radius * B;
}

/**
 * Mutual impedance per unit length between two parallel earth-return circuits (Carson, first-term
 * approximation): Zm = ωμ0/8 + j·ωμ0/(2π)·ln(De/d), De = 658.87·sqrt(ρ/f).
 * Valid for separations d much smaller than De (at 50 Hz in seawater De ≈ 45 m, in soil of 100 Ω·m ≈ 930 m);
 * for d approaching De the full Carson series or a field solution is needed.
 * At d ≥ De the first-term expression has no meaning (its reactance turns negative): the victim lies outside
 * the equivalent earth-return loop, the coupling is taken as zero and `beyondReturnDepth` is set, so a screening
 * sweep over distance runs through instead of stopping. The result is therefore discontinuous at d = De (the
 * resistive term ωμ0/8 is dropped there); do not rely on values close to De.
 * @param {object} o
 * @param {number} o.d separation of the two conductors, m
 * @param {number} o.f frequency, Hz
 * @param {number} [o.rhoEarth=1] resistivity of the return medium, Ω·m
 * @returns {{Zm:{re:number, im:number}, De:number, beyondReturnDepth:boolean}} Zm in Ω/m (zero when d ≥ De);
 *          De in m; beyondReturnDepth true when d ≥ De
 */
export function mutualImpedance({ d, f, rhoEarth = 1 }) {
  pos(d, 'd'); pos(f, 'f'); pos(rhoEarth, 'rhoEarth');
  const De = 658.87 * Math.sqrt(rhoEarth / f);
  if (d >= De) return { Zm: cx(0, 0), De, beyondReturnDepth: true };
  const w = TWO_PI * f;
  return { Zm: cx((w * MU0) / 8, ((w * MU0) / TWO_PI) * Math.log(De / d)), De, beyondReturnDepth: false };
}

/**
 * Longitudinal voltage induced in a parallel metallic pipeline or cable by the currents of a cable circuit.
 * Induced field per unit length E = −Σ Zm(d_k)·I_k over the inducing conductors (including sheath currents
 * if significant). For a balanced three-phase circuit the earth-return terms cancel and only the differences
 * in distance remain; for a fault or a DC-pole-like single current the full Zm applies.
 * Voltage: for an insulated section earthed at one end, V = E·L; for a long uniformly coated pipeline the
 * maximum voltage is limited to about E/|γ| at a discontinuity, where γ is the propagation constant of the
 * pipeline–earth circuit (argument; typically 0.1–1 per km for well-coated pipelines).
 * Conductors farther from the victim than the earth-return depth De = 658.87·sqrt(ρ/f) contribute nothing
 * (see mutualImpedance); `beyondReturnDepth` counts them.
 * @param {object} o
 * @param {{d:number, I:number|{re:number, im:number}}[]} o.conductors distance to the victim (m) and current (A)
 * @param {number} o.f frequency, Hz @param {number} [o.rhoEarth=1] Ω·m
 * @param {number} o.length length of parallelism, m
 * @param {number} [o.screeningFactor=1] reduction by earthed screening conductors (0–1)
 * @param {number} [o.gamma] magnitude of the victim's propagation constant, 1/m (for the long-line limit)
 * @returns {{E:number, voltageInsulated:number, voltageLongLine:number, voltage:number, beyondReturnDepth:number}}
 *          V/m and V (rms); `voltage` is the smaller of the two estimates when gamma is given;
 *          beyondReturnDepth = number of conductors at d ≥ De, treated as uncoupled
 */
export function inducedVoltage({ conductors, f, rhoEarth = 1, length, screeningFactor = 1, gamma }) {
  nonNeg(length, 'length');
  if (!(screeningFactor >= 0 && screeningFactor <= 1)) throw new RangeError('screeningFactor: must be in [0, 1]');
  let e = cx(0, 0); let beyondReturnDepth = 0;
  for (const c of conductors) {
    const m = mutualImpedance({ d: c.d, f, rhoEarth });
    if (m.beyondReturnDepth) beyondReturnDepth++;
    e = cadd(e, cmul(m.Zm, toCx(c.I)));
  }
  const E = cabs(e) * screeningFactor;
  const voltageInsulated = E * length;
  const voltageLongLine = gamma > 0 ? E / gamma : Infinity;
  return { E, voltageInsulated, voltageLongLine, voltage: Math.min(voltageInsulated, voltageLongLine), beyondReturnDepth };
}

/**
 * Sheath standing voltage per unit length of single-core cables with single-point-bonded (or cross-bonded,
 * within a minor section) sheaths. The voltage induced in the sheath of cable k by all conductor currents is
 * E_k = −jω·μ0/(2π)·Σ_j I_j·ln(1/d_kj), with d_kk the mean sheath radius (flux linkage of the loop sheath–remote
 * return; the remote-return terms cancel for currents summing to zero). For a trefoil this gives
 * E = ω·μ0/(2π)·I·ln(2s/d), the X·I of IEC 60287-1-1 §2.3.1.
 * Reference: IEEE Std 575 (sheath bonding), Annex on induced voltages; first-principles flux linkage.
 * @param {object} o
 * @param {{x:number, z:number, I:number|{re:number, im:number}}[]} o.conductors positions (m) and currents (A)
 * @param {number} o.sheathMeanDiameter m @param {number} o.f frequency, Hz
 * @param {number} [o.length=1] length of the single-point-bonded section, m
 * @returns {{perMetre:number[], max:number, voltage:number}} sheath voltage gradient of each cable (V/m), the
 *          largest gradient (V/m) and the largest standing voltage for the section length (V)
 */
export function sheathStandingVoltage({ conductors, sheathMeanDiameter, f, length = 1 }) {
  pos(sheathMeanDiameter, 'sheathMeanDiameter'); pos(f, 'f'); nonNeg(length, 'length');
  const w = TWO_PI * f; const rs = sheathMeanDiameter / 2;
  const perMetre = conductors.map((ck, k) => {
    let s = cx(0, 0);
    for (let j = 0; j < conductors.length; j++) {
      const cj = conductors[j];
      const d = j === k ? rs : Math.hypot(ck.x - cj.x, ck.z - cj.z);
      if (!(d > 0)) throw new RangeError('conductors: two cables coincide');
      s = cadd(s, cscale(toCx(cj.I), Math.log(1 / d)));
    }
    return cabs(s) * ((w * MU0) / TWO_PI);
  });
  const max = Math.max(...perMetre);
  return { perMetre, max, voltage: max * length };
}

/**
 * Sheath voltage limiter (SVL) selection hook: the rated voltage must exceed the power-frequency sheath voltage
 * during the worst through-fault for the fault-clearing time, U_svl ≥ k·E_fault·L.
 * The temporary-overvoltage factor k and fault gradient are project data.
 * @param {object} o
 * @param {number} o.faultGradient sheath voltage gradient under fault current, V/m (sheathStandingVoltage with the
 *        fault currents) @param {number} o.length section length, m
 * @param {number} [o.margin=1.0] ratio required between SVL rated voltage and the fault standing voltage
 * @param {number[]} [o.ratings] available SVL rated voltages, V (ascending)
 * @returns {{required:number, selected:?number}} required rated voltage (V) and the smallest adequate rating
 */
export function svlSelection({ faultGradient, length, margin = 1.0, ratings }) {
  nonNeg(faultGradient, 'faultGradient'); nonNeg(length, 'length'); pos(margin, 'margin');
  const required = faultGradient * length * margin;
  let selected = null;
  if (ratings) for (const r of ratings) if (r >= required) { selected = r; break; }
  return { required, selected };
}

/**
 * AC corrosion likelihood at a coating defect of a pipeline: current density through a circular holiday of
 * diameter dh in soil/seabed of resistivity ρ, J = 8·Vac/(ρ·π·dh) (spreading resistance of a disc, R = ρ/(2·dh)).
 * Criterion per EN 15280 / ISO 18086: corrosion unlikely below about 30 A/m² (rms) on a 1 cm² defect
 * (dh = 11.3 mm); the limit is an argument.
 * @param {object} o
 * @param {number} o.Vac AC pipe-to-soil voltage (rms), V @param {number} o.rhoSoil local resistivity, Ω·m
 * @param {number} [o.holidayDiameter=0.01128] m (1 cm² circular defect)
 * @param {number} [o.limit=30] current-density limit, A/m²
 * @returns {{J:number, limit:number, utilisation:number, ok:boolean, VacAllowable:number}} A/m², –, V
 */
export function acCorrosionCurrentDensity({ Vac, rhoSoil, holidayDiameter = 0.01128, limit = 30 }) {
  nonNeg(Vac, 'Vac'); pos(rhoSoil, 'rhoSoil'); pos(holidayDiameter, 'holidayDiameter'); pos(limit, 'limit');
  const J = (8 * Vac) / (rhoSoil * Math.PI * holidayDiameter);
  return { J, limit, utilisation: J / limit, ok: J <= limit, VacAllowable: (limit * rhoSoil * Math.PI * holidayDiameter) / 8 };
}

/**
 * Circulating current driven by an induced longitudinal voltage in a closed loop (sheath bonded at both ends,
 * pipeline earthed at both ends): I = V/|Zloop|.
 * @param {object} o
 * @param {number} o.voltage induced loop voltage, V
 * @param {{re:number, im:number}|number} o.Zloop loop impedance (self-impedance with earth return × length +
 *        earthing resistances), Ω
 * @returns {number} current magnitude, A
 */
export function circulatingCurrent({ voltage, Zloop }) {
  nonNeg(voltage, 'voltage');
  const z = cabs(toCx(Zloop));
  if (!(z > 0)) throw new RangeError('Zloop: must be non-zero');
  return voltage / z;
}
