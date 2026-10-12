/**
 * Riser stress and limit-state checks: API RP 2RD von Mises working-stress check, DNV-OS-F201 / DNV-ST-F201
 * LRFD combined-loading criteria, burst, collapse, propagating buckling, strain-based local buckling,
 * flex-joint rotation, tapered stress-joint sizing, connector capacity envelopes and fatigue hot-spot stress.
 *
 * All pressures in Pa, forces in N, moments in N·m, lengths in m, stresses in Pa. `T` in the API functions is the
 * true wall tension; `Te` in the DNV functions is the effective tension.
 */

function num(v, name, { min = -Infinity, max = Infinity, strictMin = false } = {}) {
  if (typeof v !== 'number' || !Number.isFinite(v) || v < min || v > max || (strictMin && v <= min)) {
    throw new RangeError(`risercheck: argument ${name} = ${v} is outside its valid range`);
  }
  return v;
}
function geom(D, t) {
  num(D, 'D', { min: 0, strictMin: true }); num(t, 't', { min: 0, strictMin: true, max: D / 2 });
  const Di = D - 2 * t;
  return { Di, A: Math.PI * (D * D - Di * Di) / 4, I: Math.PI * (D ** 4 - Di ** 4) / 64 };
}

/**
 * API RP 2RD von Mises check of a plain pipe under pressure, tension and bending.
 *
 * Principal stresses per API RP 2RD (1998) §5.2.3: σr = −(pe·Do + pi·Di)/(Do + Di), σθ = (pi − pe)·Do/(2t) − pi,
 * σz = T/A ± M·(Do − t)/(2I) (T = true wall tension); von Mises equivalent
 * σe = sqrt(½[(σr − σθ)² + (σθ − σz)² + (σz − σr)²]) ≤ Cf·Ca·σy with Ca = 2/3 and design case factor
 * Cf = 1.0 (operating), 1.2 (extreme and temporary), 1.35 (hydrotest), 1.5 (survival) (API RP 2RD Table 2).
 * Validity: plain pipe away from discontinuities; secondary/peak stresses not covered.
 *
 * @param {object} o
 * @param {number} o.T true wall tension [N] @param {number} [o.M=0] bending moment [N·m]
 * @param {number} [o.pi=0] internal pressure [Pa] @param {number} [o.pe=0] external pressure [Pa]
 * @param {number} o.D outer diameter [m] @param {number} o.t wall thickness [m]
 * @param {number} o.sigmaY specified minimum yield strength [Pa]
 * @param {number} [o.Cf=1.0] design case factor [-] @param {number} [o.Ca=2/3] allowable stress factor [-]
 * @returns {{ sigmaR: number, sigmaTheta: number, sigmaZ: number, vonMises: number, allowable: number,
 *   uc: number }} stresses [Pa] (sigmaZ and vonMises at the governing side of bending); uc = demand/capacity [-]
 */
export function vonMisesApi2rd({ T, M = 0, pi = 0, pe = 0, D, t, sigmaY, Cf = 1.0, Ca = 2 / 3 }) {
  const { Di, A, I } = geom(D, t);
  num(T, 'T'); num(M, 'M'); num(pi, 'pi', { min: 0 }); num(pe, 'pe', { min: 0 });
  num(sigmaY, 'sigmaY', { min: 0, strictMin: true }); num(Cf, 'Cf', { min: 0, strictMin: true }); num(Ca, 'Ca', { min: 0, strictMin: true });
  const sr = -(pe * D + pi * Di) / (D + Di);
  const sth = (pi - pe) * D / (2 * t) - pi;
  const sb = Math.abs(M) * (D - t) / (2 * I);
  const vm = (sz) => Math.sqrt(0.5 * ((sr - sth) ** 2 + (sth - sz) ** 2 + (sz - sr) ** 2));
  const z1 = T / A + sb, z2 = T / A - sb;
  const sz = vm(z1) >= vm(z2) ? z1 : z2;
  const allowable = Cf * Ca * sigmaY;
  return { sigmaR: sr, sigmaTheta: sth, sigmaZ: sz, vonMises: vm(sz), allowable, uc: vm(sz) / allowable };
}

/**
 * Burst (pressure containment) resistance, DNV-OS-F201 (2010) Eq. 5.15 / DNV-ST-F101 Eq. 5.8:
 * pb(t) = (2/√3)·(2t/(D − t))·min(fy, fu/1.15).
 * @param {object} o `{ D [m], t [m], fy [Pa], fu [Pa] }`
 * @returns {number} burst resistance [Pa]
 */
export function burstPressureDnv({ D, t, fy, fu }) {
  geom(D, t); num(fy, 'fy', { min: 0, strictMin: true }); num(fu, 'fu', { min: 0, strictMin: true });
  return (2 / Math.sqrt(3)) * (2 * t / (D - t)) * Math.min(fy, fu / 1.15);
}

/**
 * Burst check (net internal overpressure), DNV-OS-F201 (2010) Eq. 5.14: (pli − pe) ≤ pb(t1)/(γm·γSC).
 * Safety-class factor γSC = 1.04 / 1.14 / 1.26 (low / normal / high) and material factor γm = 1.15 (ULS, ALS)
 * per DNV-OS-F201 Tables 5-2 and 5-3 are the documented defaults — confirm against the governing edition.
 * @param {object} o `{ pli [Pa local incidental internal pressure], pe [Pa], D, t [m, t1 = tnom − tfab − tcorr], fy, fu
 *   [Pa], gammaSC? = 1.14, gammaM? = 1.15 }`
 * @returns {{ pb: number, demand: number, capacity: number, uc: number }} [Pa, Pa, Pa, -]
 */
export function burstCheckDnv({ pli, pe = 0, D, t, fy, fu, gammaSC = 1.14, gammaM = 1.15 }) {
  num(pli, 'pli', { min: 0 }); num(pe, 'pe', { min: 0 });
  const pb = burstPressureDnv({ D, t, fy, fu });
  const capacity = pb / (gammaM * gammaSC);
  return { pb, demand: pli - pe, capacity, uc: (pli - pe) / capacity };
}

/**
 * Burst pressure of pipe per API RP 1111 (4th ed.) Eq. 2: Pb = 0.45·(S + U)·ln(D/Di) (thick-wall form), or
 * Pb = 0.90·(S + U)·t/(D − t) (thin-wall form, D/t > 15).
 * @param {object} o `{ D [m], t [m], S [Pa SMYS], U [Pa SMTS], thin? = false }`
 * @returns {number} specified minimum burst pressure [Pa]
 */
export function burstPressureApi1111({ D, t, S, U, thin = false }) {
  const { Di } = geom(D, t); num(S, 'S', { min: 0, strictMin: true }); num(U, 'U', { min: S });
  return thin ? 0.90 * (S + U) * t / (D - t) : 0.45 * (S + U) * Math.log(D / Di);
}

/**
 * Hoop stress: thin-wall Barlow form σh = (pi − pe)·D/(2t) and mean-diameter form (pi − pe)(D − t)/(2t)
 * (ASME B31.8 §841 / DNV-ST-F101 Eq. 5.7 respectively).
 * @param {object} o `{ pi [Pa], pe [Pa], D [m], t [m] }`
 * @returns {{ barlow: number, meanDiameter: number }} [Pa]
 */
export function hoopStress({ pi, pe = 0, D, t }) {
  geom(D, t); num(pi, 'pi'); num(pe, 'pe');
  return { barlow: (pi - pe) * D / (2 * t), meanDiameter: (pi - pe) * (D - t) / (2 * t) };
}

/**
 * Characteristic collapse pressure for external overpressure, DNV-OS-F201 (2010) Eq. 5.18–5.21 /
 * DNV-ST-F101 Eq. 5.10–5.13: (pc − pel)(pc² − pp²) = pc·pel·pp·f0·D/t with pel = 2E(t/D)³/(1 − ν²),
 * pp = 2·fy·αfab·t/D and ovality f0 = (Dmax − Dmin)/D (not to be taken below 0.005). The cubic is solved with the
 * closed form of DNV-ST-F101 Sec. 13.
 * @param {object} o `{ D [m], t [m], fy [Pa], E? = 207e9 [Pa], nu? = 0.3, f0? = 0.005 [-], alphaFab? = 1.0 [-]
 *   (1.0 seamless, 0.93 UO/TRB/ERW, 0.85 UOE per DNV-ST-F101 Table 5-5) }`
 * @returns {{ pc: number, pel: number, pp: number }} collapse, elastic collapse and plastic collapse pressure [Pa]
 */
export function collapsePressureDnv({ D, t, fy, E = 207e9, nu = 0.3, f0 = 0.005, alphaFab = 1.0 }) {
  geom(D, t); num(fy, 'fy', { min: 0, strictMin: true }); num(E, 'E', { min: 0, strictMin: true });
  num(nu, 'nu', { min: 0, max: 0.499 }); num(f0, 'f0', { min: 0.005, max: 0.1 }); num(alphaFab, 'alphaFab', { min: 0.5, max: 1 });
  const pel = 2 * E * (t / D) ** 3 / (1 - nu * nu);
  const pp = 2 * fy * alphaFab * t / D;
  const b = -pel, c = -(pp * pp + pp * pel * f0 * D / t), d = pel * pp * pp;
  const u = (-b * b / 3 + c) / 3, v = 0.5 * (2 * b ** 3 / 27 - b * c / 3 + d);
  const phi = Math.acos(Math.max(-1, Math.min(1, -v / Math.sqrt(-u * u * u))));
  const y = -2 * Math.sqrt(-u) * Math.cos(phi / 3 + Math.PI / 3);
  return { pc: y - b / 3, pel, pp };
}

/**
 * Collapse check, DNV-OS-F201 (2010) Eq. 5.17: (pe − pmin) ≤ pc(t1)/(γSC·γm).
 * @param {object} o collapsePressureDnv arguments plus `{ pe [Pa], pmin? = 0 [Pa], gammaSC? = 1.14, gammaM? = 1.15 }`
 * @returns {{ pc: number, demand: number, capacity: number, uc: number }} [Pa, Pa, Pa, -]
 */
export function collapseCheckDnv(o) {
  const { pe, pmin = 0, gammaSC = 1.14, gammaM = 1.15 } = o;
  num(pe, 'pe', { min: 0 }); num(pmin, 'pmin', { min: 0 });
  const { pc } = collapsePressureDnv(o);
  const capacity = pc / (gammaSC * gammaM);
  return { pc, demand: pe - pmin, capacity, uc: (pe - pmin) / capacity };
}

/**
 * Collapse pressure per API RP 1111 Eq. 6 / API RP 2RD §5.4: Pc = Py·Pe/sqrt(Py² + Pe²), Py = 2·S·(t/D),
 * Pe = 2E·(t/D)³/(1 − ν²).
 * @param {object} o `{ D [m], t [m], S [Pa], E? = 207e9 [Pa], nu? = 0.3 }`
 * @returns {{ Pc: number, Py: number, Pe: number }} [Pa]
 */
export function collapsePressureApi({ D, t, S, E = 207e9, nu = 0.3 }) {
  geom(D, t); num(S, 'S', { min: 0, strictMin: true });
  const Py = 2 * S * t / D, Pe = 2 * E * (t / D) ** 3 / (1 - nu * nu);
  return { Pc: Py * Pe / Math.hypot(Py, Pe), Py, Pe };
}

/**
 * Propagating buckle pressure: DNV-OS-F201 (2010) Eq. 5.23 / DNV-ST-F101 Eq. 5.16 ppr = 35·fy·αfab·(t2/D)^2.5,
 * and API RP 1111 Eq. 9 Pp = 24·S·(t/D)^2.4. DNV check (Eq. 5.22): (pe − pmin) ≤ ppr/(γc·γm·γSC) with γc = 1.0
 * (no propagation allowed) or 0.9 (short-distance propagation accepted).
 * @param {object} o `{ D [m], t [m], fy [Pa], alphaFab? = 1, pe? [Pa], pmin? = 0 [Pa], gammaC? = 1.0, gammaM? = 1.15,
 *   gammaSC? = 1.14 }`
 * @returns {{ pprDnv: number, ppApi: number, capacity: number, uc: number }} [Pa, Pa, Pa, -] (uc = NaN without pe)
 */
export function propagatingBuckling({ D, t, fy, alphaFab = 1, pe, pmin = 0, gammaC = 1.0, gammaM = 1.15, gammaSC = 1.14 }) {
  geom(D, t); num(fy, 'fy', { min: 0, strictMin: true });
  const pprDnv = 35 * fy * alphaFab * (t / D) ** 2.5, ppApi = 24 * fy * (t / D) ** 2.4;
  const capacity = pprDnv / (gammaC * gammaM * gammaSC);
  return { pprDnv, ppApi, capacity, uc: pe === undefined ? NaN : (num(pe, 'pe', { min: 0 }) - pmin) / capacity };
}

/**
 * DNV-OS-F201 (2010) §5 D LRFD combined-loading criteria for pipe members under bending moment, effective
 * tension and net pressure (also DNV-ST-F201 §5.4).
 *
 * Design load effects: Md = γF·MF + γE·ME + γA·MA, Ted = γF·TeF + γE·TeE + γA·TeA (Eq. 5.5–5.6).
 * Net internal overpressure (Eq. 5.28):
 *   γSC·γm·{ (|Md|/Mk)·sqrt(1 − ((pld − pe)/pb)²) + (Ted/Tk)² } + ((pld − pe)/pb)² ≤ 1
 * Net external overpressure (Eq. 5.30):
 *   (γSC·γm)²·{ |Md|/Mk + (Ted/Tk)² }² + (γSC·γm)²·((pe − pmin)/pc)² ≤ 1
 * with plastic capacities Mk = fy·αc·(D − t2)²·t2, Tk = fy·αc·π·(D − t2)·t2, strain-hardening parameter
 * αc = (1 − β) + β·fu/fy, β = (0.4 + qh) for D/t2 < 15, (0.4 + qh)(60 − D/t2)/45 for 15 ≤ D/t2 ≤ 60, 0 above;
 * qh = ((pld − pe)/pb)·2/√3 for pld > pe, else 0; pb from burstPressureDnv(t2), pc from collapsePressureDnv(t2).
 * Default factors (ULS): γF = 1.1, γE = 1.3, γA = 1.0 (Table 5-1), γSC = 1.14 normal safety class (Table 5-2),
 * γm = 1.15 (Table 5-3) — documented defaults, to be confirmed against the governing edition.
 * Validity: D/t2 ≤ 60, ovality ≤ 3 %, |Ted| ≤ Tk; torsion and shear not included.
 * When the net internal pressure reaches the burst resistance ((pld − pe)/pb ≥ 1) the square root in Eq. 5.28 has
 * no real value: the pipe has no bending capacity left. The utilisation is then continued as
 * γSC·γm·(Ted/Tk)² + ((pld − pe)/pb)², the limit of Eq. 5.28 as the pressure ratio reaches 1 (finite, ≥ 1, growing
 * with the overpressure) and `burstExceeded` is true — the member fails on pressure containment
 * (see burstPressureDnv / the burst check), whatever the bending.
 *
 * @param {object} o
 * @param {number} o.D @param {number} o.t outer diameter and wall thickness t2 = tnom − tcorr [m]
 * @param {number} o.fy @param {number} o.fu design yield and tensile strength [Pa]
 * @param {number} [o.MF=0] @param {number} [o.ME=0] @param {number} [o.MA=0] functional / environmental /
 *   accidental bending moments [N·m]
 * @param {number} [o.TeF=0] @param {number} [o.TeE=0] @param {number} [o.TeA=0] effective tensions [N]
 * @param {number} [o.pld=0] local internal design pressure [Pa] @param {number} [o.pe=0] external pressure [Pa]
 * @param {number} [o.pmin=0] local minimum internal pressure [Pa]
 * @param {number} [o.gammaF=1.1] @param {number} [o.gammaE=1.3] @param {number} [o.gammaA=1.0]
 * @param {number} [o.gammaSC=1.14] @param {number} [o.gammaM=1.15]
 * @param {number} [o.E=207e9] @param {number} [o.nu=0.3] @param {number} [o.f0=0.005] @param {number} [o.alphaFab=1]
 * @returns {{ Md: number, Ted: number, Mk: number, Tk: number, alphaC: number, pb: number, pc: number,
 *   mode: 'internal' | 'external', uc: number, burstExceeded: boolean }} Md, Mk [N·m]; Ted, Tk [N]; pb, pc [Pa];
 *   uc is the left-hand side of the governing interaction equation (≤ 1 passes; always finite) [-]; burstExceeded =
 *   net internal pressure at or above the burst resistance
 */
export function combinedLoadingDnv(o) {
  const {
    D, t, fy, fu, MF = 0, ME = 0, MA = 0, TeF = 0, TeE = 0, TeA = 0, pld = 0, pe = 0, pmin = 0,
    gammaF = 1.1, gammaE = 1.3, gammaA = 1.0, gammaSC = 1.14, gammaM = 1.15, E = 207e9, nu = 0.3, f0 = 0.005, alphaFab = 1,
  } = o ?? {};
  geom(D, t);
  num(fy, 'fy', { min: 0, strictMin: true }); num(fu, 'fu', { min: fy });
  num(pld, 'pld', { min: 0 }); num(pe, 'pe', { min: 0 }); num(pmin, 'pmin', { min: 0 });
  if (D / t > 60) throw new RangeError('risercheck: argument D/t exceeds 60, outside the DNV combined-loading criterion');
  const Md = gammaF * MF + gammaE * ME + gammaA * MA;
  const Ted = gammaF * TeF + gammaE * TeE + gammaA * TeA;
  const pb = burstPressureDnv({ D, t, fy, fu });
  const qh = pld > pe ? ((pld - pe) / pb) * 2 / Math.sqrt(3) : 0;
  const dt = D / t;
  const beta = dt < 15 ? 0.4 + qh : (0.4 + qh) * (60 - dt) / 45;
  const alphaC = (1 - beta) + beta * fu / fy;
  const Mk = fy * alphaC * (D - t) ** 2 * t, Tk = fy * alphaC * Math.PI * (D - t) * t;
  const { pc } = collapsePressureDnv({ D, t, fy, E, nu, f0, alphaFab });
  const g = gammaSC * gammaM;
  let uc, mode, burstExceeded = false;
  if (pld >= pe) {
    mode = 'internal';
    const pr = (pld - pe) / pb;
    if (pr >= 1) { burstExceeded = true; uc = g * (Ted / Tk) ** 2 + pr * pr; }
    else uc = g * ((Math.abs(Md) / Mk) * Math.sqrt(1 - pr * pr) + (Ted / Tk) ** 2) + pr * pr;
  } else {
    mode = 'external';
    uc = g * g * (Math.abs(Md) / Mk + (Ted / Tk) ** 2) ** 2 + g * g * ((pe - pmin) / pc) ** 2;
  }
  return { Md, Ted, Mk, Tk, alphaC, pb, pc, mode, uc, burstExceeded };
}

/**
 * Strain-based local buckling (displacement-controlled bending), DNV-ST-F101 Eq. 5.30–5.31:
 * εd ≤ εc/γε, εc = 0.78·(t/D − 0.01)·(1 + 5.75·(pmin − pe)/pb)·αh^(−1.5)·αgw for pmin ≥ pe (pressure term taken as
 * zero under net external pressure, where the combined external-pressure criterion applies instead).
 * Strain resistance factor γε = 2.0 / 2.5 / 3.3 (low / medium / high safety class, DNV-ST-F101 Table 5-…) is the
 * documented default set — pass the value for the project. αh = (Rt0.5/Rm) maximum yield-to-tensile ratio
 * (typically ≤ 0.93), αgw = girth-weld factor (1.0 for D/t ≤ 20).
 * Validity: D/t ≤ 45, pi ≥ pe.
 * @param {object} o `{ strain [-] design compressive bending strain, D, t [m], alphaH? = 0.93, alphaGw? = 1,
 *   pmin? = 0, pe? = 0, pb? [Pa, needed when pmin > pe], gammaEps? = 2.5 }`
 * @returns {{ epsC: number, capacity: number, uc: number }} characteristic and design strain capacity [-], uc [-]
 */
export function localBucklingStrain({ strain, D, t, alphaH = 0.93, alphaGw = 1, pmin = 0, pe = 0, pb, gammaEps = 2.5 }) {
  geom(D, t); num(strain, 'strain', { min: 0 });
  num(alphaH, 'alphaH', { min: 0.5, max: 1 }); num(alphaGw, 'alphaGw', { min: 0, strictMin: true, max: 1 });
  num(gammaEps, 'gammaEps', { min: 1 });
  if (D / t > 45) throw new RangeError('risercheck: argument D/t exceeds 45, outside the strain criterion');
  if (t / D <= 0.01) throw new RangeError('risercheck: argument t/D must exceed 0.01');
  let pterm = 1;
  if (pmin > pe) pterm = 1 + 5.75 * (pmin - pe) / num(pb, 'pb', { min: 0, strictMin: true });
  const epsC = 0.78 * (t / D - 0.01) * pterm * alphaH ** -1.5 * alphaGw;
  return { epsC, capacity: epsC / gammaEps, uc: strain / (epsC / gammaEps) };
}

/**
 * Flex-joint rotation utilisation and restoring moment: uc = |θ|/θallow, M = kRot·θ (linear elastomeric joint).
 * @param {object} o `{ angleDeg [deg], allowableDeg [deg], kRot? [N·m/rad] }`
 * @returns {{ uc: number, moment: number }} [-], [N·m] (0 without kRot)
 */
export function flexJointCheck({ angleDeg, allowableDeg, kRot = 0 }) {
  num(angleDeg, 'angleDeg'); num(allowableDeg, 'allowableDeg', { min: 0, strictMin: true }); num(kRot, 'kRot', { min: 0 });
  return { uc: Math.abs(angleDeg) / allowableDeg, moment: kRot * Math.abs(angleDeg) * Math.PI / 180 };
}

/**
 * Tapered stress-joint sizing by the constant-curvature method (first-principles tensioned-beam statics).
 *
 * A joint clamped at its root and rotated by Δθ relative to the cable direction of the riser is tapered so that
 * its curvature κ0 is uniform. With ξ measured from the tip (junction with the plain pipe) towards the root and
 * λp = sqrt(EIp/T): the pipe beyond the tip takes the rotation κ0·λp (boundary layer M = Mtip·e^(−s/λp)), and
 * M'' = T·κ0 along the joint gives the required stiffness
 *   EI(ξ) = EIp·(1 + ξ/λp + ½(ξ/λp)²),   length L = Δθ/κ0 − λp,
 * with κ0 set by the allowable bending stress at the root outer fibre, κ0 = σb,allow/(E·rRoot) (iterated, bore
 * constant). Validity: small rotations, tension constant over the joint, no lateral load on the joint; the
 * resulting profile should be verified with the rod model.
 *
 * @param {object} o
 * @param {number} o.deltaAngle rotation to accommodate [rad] @param {number} o.T effective tension at the joint [N]
 * @param {number} o.D pipe outer diameter [m] @param {number} o.t pipe wall thickness [m]
 * @param {number} [o.E=207e9] [Pa] @param {number} o.sigmaBendAllow allowable bending stress in the joint [Pa]
 * @param {number} [o.nPoints=11]
 * @returns {{ length: number, curvature: number, lambdaPipe: number, rootOD: number, rootEI: number,
 *   rootMoment: number, xi: number[], EI: number[], OD: number[] }} length, lambdaPipe, rootOD [m]; curvature
 *   [1/m]; rootEI [N·m²]; rootMoment [N·m]; profile from tip (ξ = 0) to root [m, N·m², m]. length = 0 means the
 *   plain pipe already accommodates the rotation.
 */
export function stressJointTaper({ deltaAngle, T, D, t, E = 207e9, sigmaBendAllow, nPoints = 11 }) {
  const { Di, I } = geom(D, t);
  num(deltaAngle, 'deltaAngle', { min: 0, strictMin: true }); num(T, 'T', { min: 0, strictMin: true });
  num(sigmaBendAllow, 'sigmaBendAllow', { min: 0, strictMin: true });
  const EIp = E * I, lam = Math.sqrt(EIp / T);
  const odFromEI = (EIv) => (64 * EIv / (E * Math.PI) + Di ** 4) ** 0.25;
  let Dr = D, L = 0, k0 = 0;
  for (let it = 0; it < 200; it++) {
    k0 = sigmaBendAllow / (E * Dr / 2);
    L = Math.max(0, deltaAngle / k0 - lam);
    const q = L / lam;
    const Dn = odFromEI(EIp * (1 + q + 0.5 * q * q));
    if (Math.abs(Dn - Dr) < 1e-12 * D) { Dr = Dn; break; }
    Dr = 0.5 * (Dr + Dn);
  }
  const xi = [], EIa = [], OD = [];
  for (let k = 0; k < nPoints; k++) {
    const x = L * k / (nPoints - 1), q = x / lam;
    xi.push(x); EIa.push(EIp * (1 + q + 0.5 * q * q)); OD.push(odFromEI(EIa[k]));
  }
  const rootEI = EIa[nPoints - 1];
  return { length: L, curvature: k0, lambdaPipe: lam, rootOD: OD[nPoints - 1], rootEI, rootMoment: rootEI * k0, xi, EI: EIa, OD };
}

/**
 * Connector / hang-off load utilisation against a rated capacity envelope given as a convex polygon in the
 * (tension, bending moment) plane (e.g. a manufacturer capacity chart at the relevant pressure). The utilisation
 * is the ratio of the load vector to the envelope along the same ray from the origin.
 * @param {object} o `{ T [N], M [N·m], envelope: Array<{ T [N], M [N·m] }> }` vertices in order around the origin
 *   (at least 3; the origin must lie inside)
 * @returns {{ uc: number, capacityT: number, capacityM: number }} [-]; envelope point on the load ray [N, N·m]
 */
export function connectorUtilisation({ T, M, envelope }) {
  num(T, 'T'); num(M, 'M');
  if (!Array.isArray(envelope) || envelope.length < 3) throw new RangeError('risercheck: argument envelope needs at least 3 vertices');
  if (T === 0 && M === 0) return { uc: 0, capacityT: 0, capacityM: 0 };
  let best = Infinity;
  for (let k = 0; k < envelope.length; k++) {
    const a = envelope[k], b = envelope[(k + 1) % envelope.length];
    // solve s·(T, M) = a + u·(b − a)
    const ex = b.T - a.T, ey = b.M - a.M, det = -T * ey + M * ex;
    if (Math.abs(det) < 1e-300) continue;
    const s = (-a.T * ey + a.M * ex) / det, u = (T * a.M - M * a.T) / det;
    if (s > 0 && u >= -1e-12 && u <= 1 + 1e-12 && s < best) best = s;
  }
  if (!Number.isFinite(best)) throw new RangeError('risercheck: argument envelope does not enclose the origin');
  return { uc: 1 / best, capacityT: best * T, capacityM: best * M };
}

/**
 * Fatigue hot-spot stress range from tension and bending-moment ranges and the resulting Miner damage with a
 * caller-supplied S–N function: Δσ = SCF·(ΔT/A + ΔM·r/I), damage = n/N(Δσ) (DNV-RP-C203 §2; DNV-OS-F201 §5 F).
 * @param {object} o `{ dT [N], dM [N·m], D [m], t [m], scf? = 1, radius? [m, default D/2 outer fibre], cycles? = 1,
 *   snCyclesToFailure?: (stressRangePa) => N }`
 * @returns {{ stressRange: number, axialRange: number, bendingRange: number, damage: number }} [Pa]; damage [-]
 *   (NaN when no S–N function is given)
 */
export function fatigueHotSpot({ dT = 0, dM = 0, D, t, scf = 1, radius, cycles = 1, snCyclesToFailure }) {
  const { A, I } = geom(D, t);
  num(dT, 'dT'); num(dM, 'dM'); num(scf, 'scf', { min: 0, strictMin: true }); num(cycles, 'cycles', { min: 0 });
  const r = radius ?? D / 2;
  const ax = Math.abs(dT) / A, bd = Math.abs(dM) * r / I;
  const S = scf * (ax + bd);
  let damage = NaN;
  if (typeof snCyclesToFailure === 'function') {
    const N = snCyclesToFailure(S);
    damage = N > 0 && Number.isFinite(N) ? cycles / N : (N === Infinity ? 0 : NaN);
  }
  return { stressRange: S, axialRange: ax, bendingRange: bd, damage };
}
