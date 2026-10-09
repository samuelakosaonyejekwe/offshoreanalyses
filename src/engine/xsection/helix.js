/**
 * Mechanics of helically wound layers (armour wires, umbilical components, cable armour).
 *
 * Conventions used throughout this file
 * -------------------------------------
 *  - alpha  lay angle measured from the pipe axis, rad, 0 < alpha < pi/2.
 *  - R      helix (pitch-circle) radius, m.
 *  - lay    +1 or -1: hand of the helix. A positive twist per unit length elongates wires with lay = +1.
 *  - psi    angular position of a wire on the cross-section, rad, measured from the extrados (the fibre with the
 *           largest tensile bending strain), so the bending strain of the supporting cylinder is kappa*R*cos(psi).
 *  - A wire has axial stiffness E*A, bending stiffness E*In about the axis lying in the cylinder surface (normal
 *    curvature, i.e. bending "over" the cylinder, governed by the radial wire thickness), E*It about the surface
 *    normal (transverse curvature, governed by the wire width) and torsional stiffness G*J.
 *
 * Kinematics (first principles). A wire element of length dl has dz = dl*cos(alpha), R*dtheta = dl*sin(alpha).
 * Stretching the cylinder by eps, twisting it by tau (rad/m) and expanding the radius by dR gives
 *     eps_w = cos^2(alpha)*eps + lay*R*sin(alpha)*cos(alpha)*tau + sin^2(alpha)*dR/R              (K1)
 * which is the small-strain projection used by Feret & Bournazel (1987), Witz & Tan (1992) and Saevik (1992).
 */

const TWO_PI = 2 * Math.PI;

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
function angle(v, name = 'alpha') {
  if (!(typeof v === 'number' && v > 0 && v < Math.PI / 2)) throw new RangeError(`${name} must lie strictly between 0 and pi/2 rad`);
  return v;
}
function laySign(v) {
  if (v === undefined) return 1;
  if (v !== 1 && v !== -1) throw new RangeError('lay must be +1 or -1');
  return v;
}

/**
 * Geometry of a helical layer.
 * Pitch L = 2*pi*R/tan(alpha); wire length per pitch L/cos(alpha); a wire of width b (measured normal to its own
 * axis in the cylinder surface) occupies b/cos(alpha) of the circumference, hence the fill fraction
 * F = n*b/(2*pi*R*cos(alpha)). The unbent helix has curvature sin^2(alpha)/R and torsion sin(alpha)cos(alpha)/R.
 * Reference: differential geometry of the circular helix; API RP 17B §5 terminology.
 * Validity: 0 < alpha < pi/2. Give either nWires or fillFraction (0 < F <= 1); with fillFraction the largest
 * whole number of wires not exceeding it is returned.
 * @param {object} o
 * @param {number} o.R helix radius, m
 * @param {number} o.alpha lay angle, rad
 * @param {number} o.wireWidth wire width b in the cylinder surface, m
 * @param {number} [o.nWires] number of wires, -
 * @param {number} [o.fillFraction] target fill fraction, -
 * @returns {{pitch:number, wireLengthPerPitch:number, wireLengthPerMetre:number, nWires:number,
 *   fillFraction:number, gap:number, curvature:number, torsion:number}} pitch m; wire lengths m and m/m;
 *   gap = clear gap between adjacent wires measured normal to the wire axis, m; curvature 1/m; torsion rad/m
 */
export function helixGeometry({ R, alpha, wireWidth, nWires, fillFraction }) {
  pos(R, 'R'); angle(alpha); pos(wireWidth, 'wireWidth');
  const c = Math.cos(alpha), s = Math.sin(alpha);
  const maxN = TWO_PI * R * c / wireWidth;
  let n;
  if (nWires !== undefined) {
    n = pos(nWires, 'nWires');
    if (n > maxN * (1 + 1e-12)) throw new RangeError('nWires exceeds the number that fits on the circumference');
  } else {
    if (!(fillFraction > 0 && fillFraction <= 1)) throw new RangeError('fillFraction must be in (0, 1]');
    n = Math.floor(fillFraction * maxN + 1e-9);
    if (n < 1) throw new RangeError('fillFraction too small for one wire');
  }
  const pitch = TWO_PI * R / Math.tan(alpha);
  return {
    pitch,
    wireLengthPerPitch: pitch / c,
    wireLengthPerMetre: 1 / c,
    nWires: n,
    fillFraction: n / maxN,
    gap: TWO_PI * R * c / n - wireWidth,
    curvature: s * s / R,
    torsion: s * c / R,
  };
}

/**
 * Section properties of a rectangular wire (width b in the cylinder surface, radial thickness t).
 * In = b t^3/12 (normal curvature), It = t b^3/12 (transverse curvature); Saint-Venant torsion constant
 * J = a c^3 [1/3 - 0.21 (c/a)(1 - c^4/(12 a^4))] with a the longer and c the shorter side
 * (Roark's Formulas for Stress and Strain, Table 10.1, solid rectangle).
 * @param {object} o
 * @param {number} o.width b, m
 * @param {number} o.thickness t, m
 * @returns {{A:number, In:number, It:number, J:number}} m^2, m^4, m^4, m^4
 */
export function rectWireSection({ width, thickness }) {
  const b = pos(width, 'width'), t = pos(thickness, 'thickness');
  const a = Math.max(b, t), c = Math.min(b, t);
  return {
    A: b * t,
    In: b * t ** 3 / 12,
    It: t * b ** 3 / 12,
    J: a * c ** 3 * (1 / 3 - 0.21 * (c / a) * (1 - c ** 4 / (12 * a ** 4))),
  };
}

/**
 * Section properties of a round wire of diameter d: A = pi d^2/4, In = It = pi d^4/64, J = pi d^4/32.
 * @param {object} o
 * @param {number} o.diameter d, m
 * @returns {{A:number, In:number, It:number, J:number}} m^2, m^4, m^4, m^4
 */
export function roundWireSection({ diameter }) {
  const d = pos(diameter, 'diameter');
  const I = Math.PI * d ** 4 / 64;
  return { A: Math.PI * d * d / 4, In: I, It: I, J: 2 * I };
}

/**
 * Axial strain of a helical wire from the deformation of its supporting cylinder, Eq. (K1) of the file header.
 * Validity: small strains; wire follows the cylinder surface (no gap formation).
 * @param {object} o
 * @param {number} o.alpha lay angle, rad
 * @param {number} o.R helix radius, m
 * @param {number} [o.lay=1] hand of lay, +1/-1
 * @param {number} [o.epsAxial=0] cylinder axial strain, -
 * @param {number} [o.twist=0] cylinder twist per unit length, rad/m
 * @param {number} [o.dR=0] radial displacement of the helix radius, m (positive outward)
 * @returns {{strain:number}} wire axial strain, -
 */
export function wireStrain({ alpha, R, lay, epsAxial = 0, twist = 0, dR = 0 }) {
  angle(alpha); pos(R, 'R'); const sg = laySign(lay);
  fin(epsAxial, 'epsAxial'); fin(twist, 'twist'); fin(dR, 'dR');
  const c = Math.cos(alpha), s = Math.sin(alpha);
  return { strain: c * c * epsAxial + sg * R * s * c * twist + s * s * dR / R };
}

/**
 * Stiffness contributions of one helical layer of n wires, from (K1) and the resolution of the wire force
 * Fw = E A eps_w into an axial force n Fw cos(alpha), a torque lay n Fw R sin(alpha) and a hoop force per unit
 * length n Fw sin^2(alpha)/(2 pi R cos(alpha)) (n tan(alpha)/(2 pi R) wires cross a longitudinal cut per metre,
 * each carrying Fw sin(alpha)):
 *     kAxial         = n E A cos^3(alpha)                         N          (axial force per unit axial strain)
 *     kTorsion       = n E A R^2 sin^2(alpha) cos(alpha)          N m^2      (torque per unit twist rad/m)
 *     kCoupling      = lay n E A R sin(alpha) cos^2(alpha)        N m        (force per twist = torque per strain)
 *     kAxialRadial   = n E A sin^2(alpha) cos(alpha)              N          (axial force per unit dR/R)
 *     kTorsionRadial = lay n E A R sin^3(alpha)                   N m        (torque per unit dR/R)
 *     kRadial        = n E A sin^4(alpha)/cos(alpha)              N          (2 pi R x hoop force/length per dR/R)
 * The radial terms are the thin-layer contraction coupling: the generalised force conjugate to dR/R is
 * 2 pi R N_theta, so the 3x3 matrix in (eps, tau, dR/R) is symmetric.
 * Reference: Witz & Tan (1992) Marine Structures 5, 205-227; Feret & Bournazel (1987) J. Offshore Mech. Arct.
 * Eng. 109, 263-269.
 * Validity: linear elastic wires, small strain, wire bending/torsion stiffness neglected in the axisymmetric
 * response.
 * @param {object} o
 * @param {number} o.n number of wires, -
 * @param {number} o.E Young's modulus, Pa
 * @param {number} o.A wire area, m^2
 * @param {number} o.R helix radius, m
 * @param {number} o.alpha lay angle, rad
 * @param {number} [o.lay=1] hand of lay
 * @returns {{kAxial:number, kTorsion:number, kCoupling:number, kAxialRadial:number, kTorsionRadial:number,
 *   kRadial:number}} units as listed above
 */
export function helixLayerStiffness({ n, E, A, R, alpha, lay }) {
  pos(n, 'n'); pos(E, 'E'); pos(A, 'A'); pos(R, 'R'); angle(alpha); const sg = laySign(lay);
  const c = Math.cos(alpha), s = Math.sin(alpha), k = n * E * A;
  return {
    kAxial: k * c ** 3,
    kTorsion: k * R * R * s * s * c,
    kCoupling: sg * k * R * s * c * c,
    kAxialRadial: k * s * s * c,
    kTorsionRadial: sg * k * R * s ** 3,
    kRadial: k * s ** 4 / c,
  };
}

/**
 * Contact pressure exerted on the underlying layer by tensioned helical wires.
 * Derivation: a wire with tension Fw and curvature sin^2(alpha)/R needs an inward line load
 * w = Fw sin^2(alpha)/R per unit wire length. Per metre of pipe each wire has length 1/cos(alpha), so the n wires
 * give a radial force n Fw sin^2(alpha)/(R cos(alpha)) per metre, spread over the circumference 2 pi R:
 *     q = n Fw sin^2(alpha) / (2 pi R^2 cos(alpha)),  Fw = sigma A.
 * Reference: Feret & Bournazel (1987); identical to hoop equilibrium N_theta = q R.
 * Validity: layer-average (smeared) pressure; the pressure on the wire footprint is q / fillFraction.
 * @param {object} o
 * @param {number} o.n number of wires
 * @param {number} o.wireForce wire tension Fw, N (negative = compression gives outward load)
 * @param {number} o.R helix radius, m
 * @param {number} o.alpha lay angle, rad
 * @returns {{pressure:number, lineLoadPerWire:number, radialForcePerLength:number}} Pa (smeared, positive
 *   inward); N/m per unit wire length; N/m per unit pipe length
 */
export function helixContactPressure({ n, wireForce, R, alpha }) {
  pos(n, 'n'); fin(wireForce, 'wireForce'); pos(R, 'R'); angle(alpha);
  const c = Math.cos(alpha), s = Math.sin(alpha);
  const lineLoadPerWire = wireForce * s * s / R;
  const radialForcePerLength = n * lineLoadPerWire / c;
  return { pressure: radialForcePerLength / (TWO_PI * R), lineLoadPerWire, radialForcePerLength };
}

/**
 * Change of wire curvature and twist caused by bending the supporting cylinder to curvature kappa, for the two
 * limiting wire paths. Derived to first order in R*kappa on the torus obtained by bending the cylinder:
 *
 *  - surface principal curvatures: 1/R (hoop) and kappa cos(psi) (longitudinal); normal curvature of a curve at
 *    angle a to the axis (Euler) kn = sin^2(a)/R + cos^2(a) kappa cos(psi); geodesic torsion
 *    tg = (1/R - kappa cos(psi)) sin(a) cos(a); geodesic curvature (Liouville, surface of revolution)
 *    kg = -da/ds - cos(a) kappa sin(psi).
 *  - loxodromic path (wire stays on the material line of the cylinder, no slip): the lay angle changes by
 *    da = -sin(alpha)cos(alpha) R kappa cos(psi) and the wire stretches by cos^2(alpha) R kappa cos(psi). Expressed
 *    per unit undeformed wire length:
 *        dKn = cos^4(alpha) kappa cos(psi)
 *        dKt = -(1 + sin^2(alpha)) cos(alpha) kappa sin(psi)
 *        dTw = -sin(alpha) cos^3(alpha) kappa cos(psi)
 *    (the classical expressions quoted by Saevik 1992/2011 and Witz & Tan 1992).
 *  - geodesic path (wire free to slip sideways to the shortest path; Clairaut's relation
 *    (1 + R kappa cos(psi)) cos(a) = const gives da = cot(alpha) R kappa cos(psi)):
 *        dKn = 3 cos^2(alpha) kappa cos(psi)
 *        dKt = 0
 *        dTw = cot(alpha) (1 - 3 sin^2(alpha)) kappa cos(psi)
 *    (Feret & Bournazel 1987 give the same factor 3 on the normal curvature).
 * Real armour wires lie between the two; the loxodromic path is the usual choice for wires restrained sideways
 * by friction and neighbours, the geodesic an upper bound on normal-curvature bending stress.
 * Validity: R*kappa << 1; signs follow the header convention (normal curvature positive = same sense as 1/R).
 * @param {object} o
 * @param {number} o.alpha lay angle, rad
 * @param {number} o.kappa cylinder curvature, 1/m
 * @param {number} [o.psi=0] angular position from the extrados, rad
 * @param {'loxodromic'|'geodesic'} [o.curve='loxodromic'] wire path assumption
 * @returns {{dKn:number, dKt:number, dTw:number}} changes of normal curvature, transverse curvature (1/m) and
 *   twist (rad/m)
 */
export function wireCurvatureChange({ alpha, kappa, psi = 0, curve = 'loxodromic' }) {
  angle(alpha); fin(kappa, 'kappa'); fin(psi, 'psi');
  const c = Math.cos(alpha), s = Math.sin(alpha), cp = Math.cos(psi), sp = Math.sin(psi);
  if (curve === 'loxodromic') {
    return { dKn: c ** 4 * kappa * cp, dKt: -(1 + s * s) * c * kappa * sp, dTw: -s * c ** 3 * kappa * cp };
  }
  if (curve === 'geodesic') {
    return { dKn: 3 * c * c * kappa * cp, dKt: 0, dTw: (c / s) * (1 - 3 * s * s) * kappa * cp };
  }
  throw new RangeError("curve must be 'loxodromic' or 'geodesic'");
}

/**
 * Local bending stresses in a wire from the curvature changes of {@link wireCurvatureChange}:
 * sigmaNormal = E (t/2) dKn at the inner/outer wire faces, sigmaTransverse = E (b/2) dKt at the wire edges; the
 * corner stress is the sum of the magnitudes (Saevik 2011, Marine Structures 24, 97-111).
 * Validity: elastic wire, R*kappa << 1.
 * @param {object} o
 * @param {number} o.E Young's modulus, Pa
 * @param {number} o.thickness radial wire thickness t, m
 * @param {number} o.width wire width b, m
 * @param {number} o.alpha lay angle, rad
 * @param {number} o.kappa cylinder curvature, 1/m
 * @param {number} [o.psi=0] angular position, rad
 * @param {'loxodromic'|'geodesic'} [o.curve='loxodromic']
 * @returns {{sigmaNormal:number, sigmaTransverse:number, sigmaCorner:number, dKn:number, dKt:number,
 *   dTw:number}} stresses Pa (signed, corner = sum of magnitudes), curvatures 1/m
 */
export function wireBendingStress({ E, thickness, width, alpha, kappa, psi = 0, curve = 'loxodromic' }) {
  pos(E, 'E'); pos(thickness, 'thickness'); pos(width, 'width');
  const k = wireCurvatureChange({ alpha, kappa, psi, curve });
  const sigmaNormal = E * thickness / 2 * k.dKn;
  const sigmaTransverse = E * width / 2 * k.dKt;
  return { sigmaNormal, sigmaTransverse, sigmaCorner: Math.abs(sigmaNormal) + Math.abs(sigmaTransverse), ...k };
}

/**
 * Bending stiffness contribution of one helical layer in the full-stick and full-slip states.
 *
 * Full stick (plane sections, wires bonded to the cylinder): wire strain cos^2(alpha) R kappa cos(psi), force
 * component along the pipe Fw cos(alpha) at lever arm R cos(psi); summing n wires spread round the circumference
 * (mean of cos^2 = 1/2):
 *     EI_stick = n E A R^2 cos^3(alpha) / 2.
 * Written with the fill fraction F and layer thickness t (n A = F 2 pi R t cos(alpha)) this is
 * F pi R^3 t E cos^4(alpha), which is where the cos^4 form found in the literature comes from.
 *
 * Full slip (no axial force build-up; each wire bends and twists about its own axes). Strain energy per metre of
 * pipe of n wires of length 1/cos(alpha) with the curvature changes of {@link wireCurvatureChange}:
 *     loxodromic: EI_slip = n/(2 cos a) [ E In cos^8 a + E It (1+sin^2 a)^2 cos^2 a + G J sin^2 a cos^6 a ]
 *     geodesic:   EI_slip = n/(2 cos a) [ 9 E In cos^4 a + G J cot^2 a (1 - 3 sin^2 a)^2 ]
 * Reference: Witz & Tan (1992) Marine Structures 5, 229-249; Saevik (1992) thesis, NTH.
 * Validity: R*kappa << 1, wires evenly distributed (n >= 3).
 * @param {object} o
 * @param {number} o.n number of wires
 * @param {number} o.E Young's modulus, Pa
 * @param {number} o.A wire area, m^2
 * @param {number} o.R helix radius, m
 * @param {number} o.alpha lay angle, rad
 * @param {number} [o.In=0] wire second moment for normal curvature, m^4
 * @param {number} [o.It=0] wire second moment for transverse curvature, m^4
 * @param {number} [o.J=0] wire torsion constant, m^4
 * @param {number} [o.G] shear modulus, Pa (default E/2.6, i.e. nu = 0.3)
 * @param {'loxodromic'|'geodesic'} [o.curve='loxodromic']
 * @returns {{stick:number, slip:number}} N m^2
 */
export function helixBendingStiffness({ n, E, A, R, alpha, In = 0, It = 0, J = 0, G, curve = 'loxodromic' }) {
  pos(n, 'n'); pos(E, 'E'); pos(A, 'A'); pos(R, 'R'); angle(alpha);
  nonneg(In, 'In'); nonneg(It, 'It'); nonneg(J, 'J');
  const Gs = G === undefined ? E / 2.6 : pos(G, 'G');
  const c = Math.cos(alpha), s = Math.sin(alpha);
  const stick = n * E * A * R * R * c ** 3 / 2;
  let slip;
  if (curve === 'loxodromic') {
    slip = n / (2 * c) * (E * In * c ** 8 + E * It * (1 + s * s) ** 2 * c * c + Gs * J * s * s * c ** 6);
  } else if (curve === 'geodesic') {
    slip = n / (2 * c) * (9 * E * In * c ** 4 + Gs * J * (c / s) ** 2 * (1 - 3 * s * s) ** 2);
  } else throw new RangeError("curve must be 'loxodromic' or 'geodesic'");
  return { stick, slip };
}

/**
 * Slip-onset (critical) curvature of a helical wire held by friction, and the associated friction stress and
 * moment, from the shear-lag balance along the wire.
 *
 * Derivation. In stick the wire force is Q = E A cos^2(alpha) R kappa cos(psi). Along the wire
 * d(psi)/ds = sin(alpha)/R, so |dQ/ds| is largest at the neutral axis (psi = pi/2):
 * |dQ/ds|max = E A cos^2(alpha) sin(alpha) kappa. Slip starts when this equals the available friction per unit
 * wire length f:
 *     kappa_cr = f / (E A cos^2(alpha) sin(alpha)).
 * With further bending the slip zone spreads from the neutral axis towards extrados and intrados; inside it
 * Q = (f R / sin(alpha)) (pi/2 - psi). Matching the stick solution at the stick/slip boundary psi0 gives
 * kappa/kappa_cr = (pi/2 - psi0)/cos(psi0); full slip (psi0 = 0) is reached at kappa = (pi/2) kappa_cr, where the
 * wire force at the extrados is the friction force collected over a quarter pitch:
 *     sigma_f,max = pi R f / (2 sin(alpha) A).
 * The friction moment of n fully slipped wires is the first Fourier component of the triangular force
 * distribution: M_f = (2/pi) n f R^2 cot(alpha) = (4/pi) x (stick moment at kappa_cr).
 * This is the formulation of Feret & Bournazel (1987) and Saevik (2011, Marine Structures 24, 97-111).
 *
 * f may be given directly (frictionPerLength) or built as f = mu (qOuter + qInner) x tributary width, where the q
 * are layer-average contact pressures on the outer and inner faces and the tributary width of one wire is
 * 2 pi R cos(alpha)/n.
 * Validity: Coulomb friction, constant contact pressure round the circumference, R*kappa << 1.
 * @param {object} o
 * @param {number} o.E Young's modulus, Pa
 * @param {number} o.A wire area, m^2
 * @param {number} o.R helix radius, m
 * @param {number} o.alpha lay angle, rad
 * @param {number} [o.n=1] number of wires (needed for the friction moment and for the pressure form of f)
 * @param {number} [o.frictionPerLength] available friction per unit wire length f, N/m
 * @param {number} [o.mu] friction coefficient, - (used when frictionPerLength is omitted)
 * @param {number} [o.qOuter=0] smeared contact pressure on the outer face, Pa
 * @param {number} [o.qInner=0] smeared contact pressure on the inner face, Pa
 * @returns {{frictionPerLength:number, kappaCr:number, kappaFullSlip:number, sigmaStickAtOnset:number,
 *   sigmaFrictionMax:number, momentAtOnset:number, frictionMoment:number, slipAmplitudePerCurvature:number}}
 *   N/m; 1/m; 1/m; Pa; Pa; N m; N m; m per (1/m): relative wire slip at the neutral axis of a fully slipped wire
 *   per unit curvature, R^2 cos^2(alpha)/sin(alpha)
 */
export function criticalCurvature({ E, A, R, alpha, n = 1, frictionPerLength, mu, qOuter = 0, qInner = 0 }) {
  pos(E, 'E'); pos(A, 'A'); pos(R, 'R'); angle(alpha); pos(n, 'n');
  const c = Math.cos(alpha), s = Math.sin(alpha);
  let f = frictionPerLength;
  if (f === undefined) {
    nonneg(mu, 'mu'); nonneg(qOuter, 'qOuter'); nonneg(qInner, 'qInner');
    f = mu * (qOuter + qInner) * TWO_PI * R * c / n;
  } else nonneg(f, 'frictionPerLength');
  const kappaCr = f / (E * A * c * c * s);
  const momentAtOnset = n * E * A * R * R * c ** 3 / 2 * kappaCr;
  return {
    frictionPerLength: f,
    kappaCr,
    kappaFullSlip: Math.PI / 2 * kappaCr,
    sigmaStickAtOnset: E * c * c * R * kappaCr,
    sigmaFrictionMax: Math.PI * R * f / (2 * s * A),
    momentAtOnset,
    frictionMoment: 4 / Math.PI * momentAtOnset,
    slipAmplitudePerCurvature: R * R * c * c / s,
  };
}

/**
 * Axial (friction) stress in a wire at angular position psi for monotonic bending from the straight state,
 * covering stick, partial slip and full slip (see {@link criticalCurvature} for the derivation). The stick/slip
 * boundary psi0 solves kappa/kappa_cr = (pi/2 - psi0)/cos(psi0) (bisection); for |psi| < psi0 (folded into the
 * first quadrant) the wire is in stick, sigma = E cos^2(alpha) R kappa cos(psi), otherwise
 * sigma = sigma_f,max (1 - 2 psi/pi).
 * Validity: first loading from zero curvature; for cyclic bending use ranges (twice the amplitude values,
 * Masing rule) as done in flexfatigue.js.
 * @param {object} o
 * @param {number} o.E Young's modulus, Pa
 * @param {number} o.A wire area, m^2
 * @param {number} o.R helix radius, m
 * @param {number} o.alpha lay angle, rad
 * @param {number} o.kappa curvature (>= 0), 1/m
 * @param {number} o.psi angular position, rad
 * @param {number} o.frictionPerLength available friction per unit wire length, N/m
 * @returns {{sigma:number, state:'stick'|'partial-slip'|'full-slip', psi0:number, slipped:boolean}} Pa; rad
 */
export function wireFrictionStress({ E, A, R, alpha, kappa, psi, frictionPerLength }) {
  nonneg(kappa, 'kappa'); fin(psi, 'psi');
  const cc = criticalCurvature({ E, A, R, alpha, frictionPerLength });
  const c = Math.cos(alpha);
  // fold psi into [0, pi/2] using the cos symmetry/antisymmetry of the force distribution
  let p = Math.abs(Math.atan2(Math.sin(psi), Math.cos(psi)));
  let sign = 1;
  if (p > Math.PI / 2) { p = Math.PI - p; sign = -1; }
  const stickSigma = E * c * c * R * kappa * Math.cos(p);
  if (kappa <= cc.kappaCr) return { sigma: sign * stickSigma, state: 'stick', psi0: Math.PI / 2, slipped: false };
  const slipSigma = cc.sigmaFrictionMax * (1 - 2 * p / Math.PI);
  if (kappa >= cc.kappaFullSlip) return { sigma: sign * slipSigma, state: 'full-slip', psi0: 0, slipped: true };
  const ratio = kappa / cc.kappaCr;
  let lo = 0, hi = Math.PI / 2 - 1e-12;
  for (let i = 0; i < 80; i++) {
    const mid = 0.5 * (lo + hi);
    // (pi/2 - x)/cos x decreases from pi/2 at x = 0 to 1 at x = pi/2
    if ((Math.PI / 2 - mid) / Math.cos(mid) > ratio) lo = mid; else hi = mid;
  }
  const psi0 = 0.5 * (lo + hi);
  const slipped = p > psi0;
  return { sigma: sign * (slipped ? slipSigma : stickSigma), state: 'partial-slip', psi0, slipped };
}

/**
 * Total wire stress at an angular position: mean axial stress from the axisymmetric loads, friction stress
 * ({@link wireFrictionStress}) and local bending stresses ({@link wireBendingStress}).
 * @param {object} o
 * @param {number} o.E Pa
 * @param {number} o.A wire area, m^2
 * @param {number} o.R helix radius, m
 * @param {number} o.alpha lay angle, rad
 * @param {number} o.thickness radial wire thickness, m
 * @param {number} o.width wire width, m
 * @param {number} o.kappa curvature (>= 0), 1/m
 * @param {number} o.psi angular position, rad
 * @param {number} o.frictionPerLength N/m
 * @param {number} [o.sigmaAxial=0] wire stress from tension/pressure/torque, Pa
 * @param {'loxodromic'|'geodesic'} [o.curve='loxodromic']
 * @returns {{sigmaAxial:number, sigmaFriction:number, sigmaNormal:number, sigmaTransverse:number,
 *   sigmaMaxCorner:number, sigmaMinCorner:number, state:string}} Pa; extreme corner stresses of the wire section
 */
export function wireStressAtPosition({ E, A, R, alpha, thickness, width, kappa, psi, frictionPerLength, sigmaAxial = 0, curve = 'loxodromic' }) {
  fin(sigmaAxial, 'sigmaAxial');
  const fr = wireFrictionStress({ E, A, R, alpha, kappa, psi, frictionPerLength });
  const b = wireBendingStress({ E, thickness, width, alpha, kappa, psi, curve });
  const membrane = sigmaAxial + fr.sigma;
  return {
    sigmaAxial,
    sigmaFriction: fr.sigma,
    sigmaNormal: b.sigmaNormal,
    sigmaTransverse: b.sigmaTransverse,
    sigmaMaxCorner: membrane + b.sigmaCorner,
    sigmaMinCorner: membrane - b.sigmaCorner,
    state: fr.state,
  };
}

/**
 * Moment-curvature hysteresis loop of a cross-section with friction-held helical layers, bilinear idealisation.
 * The friction-carried part of the moment is elastic-perfectly-frictional: stiffness (EIstick - EIslip) up to the
 * friction moment Mf, then constant; the slip part EIslip*kappa is always elastic. Yield (knee) curvature
 * kappaY = Mf/(EIstick - EIslip); with Mf from {@link criticalCurvature} this equals (4/pi) kappa_cr, i.e. the
 * bilinear law replaces the gradual transition between kappa_cr and (pi/2) kappa_cr by a sharp knee enclosing
 * the same friction moment. For a steady cycle of amplitude kappaAmp > kappaY (Masing behaviour) the loop is a
 * parallelogram with
 *     energy per cycle  Ed = 4 Mf (kappaAmp - kappaY)                       (J per metre of pipe)
 *     peak moment       Mmax = Mf + EIslip kappaAmp
 *     equivalent viscous damping ratio  zeta = Ed / (2 pi Mmax kappaAmp)    (Ed / (4 pi Es), Es = Mmax kappaAmp / 2)
 * Reference: bilinear hysteresis, e.g. Chopra, Dynamics of Structures §3.9 (equivalent damping); friction moment
 * after Saevik (2011).
 * Validity: constant contact pressure during the cycle; rate-independent friction.
 * @param {object} o
 * @param {number} o.EIstick stick bending stiffness of the section, N m^2
 * @param {number} o.EIslip slip bending stiffness of the section, N m^2 (< EIstick)
 * @param {number} o.frictionMoment Mf, N m
 * @param {number} o.kappaAmp curvature amplitude of the cycle, 1/m
 * @param {number} [o.nPoints=41] points per branch (>= 2)
 * @returns {{kappaY:number, Mmax:number, energyPerCycle:number, dampingRatio:number, secantStiffness:number,
 *   kappa:number[], moment:number[]}} 1/m; N m; J/m; -; N m^2; closed loop arrays starting at (+kappaAmp, +Mmax)
 */
export function momentCurvatureLoop({ EIstick, EIslip, frictionMoment, kappaAmp, nPoints = 41 }) {
  pos(EIstick, 'EIstick'); nonneg(EIslip, 'EIslip'); nonneg(frictionMoment, 'frictionMoment'); pos(kappaAmp, 'kappaAmp');
  if (!(EIslip < EIstick)) throw new RangeError('EIslip must be smaller than EIstick');
  if (!(Number.isInteger(nPoints) && nPoints >= 2)) throw new RangeError('nPoints must be an integer >= 2');
  const kf = EIstick - EIslip;
  const kappaY = frictionMoment / kf;
  const slips = kappaAmp > kappaY;
  const mfPeak = slips ? frictionMoment : kf * kappaAmp;
  const Mmax = mfPeak + EIslip * kappaAmp;
  const energyPerCycle = slips ? 4 * frictionMoment * (kappaAmp - kappaY) : 0;
  const kappa = [], moment = [];
  // unloading branch from +kappaAmp to -kappaAmp, then reloading (point symmetric)
  for (let dir = -1; dir <= 1; dir += 2) {
    for (let i = 0; i < nPoints; i++) {
      const k = -dir * kappaAmp + dir * 2 * kappaAmp * i / (nPoints - 1);
      const travelled = Math.abs(k + dir * kappaAmp);
      const mf = -dir * mfPeak + dir * Math.min(kf * travelled, 2 * mfPeak);
      kappa.push(k); moment.push(mf + EIslip * k);
    }
  }
  return {
    kappaY, Mmax, energyPerCycle,
    dampingRatio: energyPerCycle / (TWO_PI * Mmax * kappaAmp),
    secantStiffness: Mmax / kappaAmp,
    kappa, moment,
  };
}
