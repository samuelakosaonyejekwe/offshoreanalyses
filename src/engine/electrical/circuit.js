/**
 * Cable and transmission circuit calculations: line constants, sequence impedances, lumped and distributed line
 * models, reactive compensation, HVDC links, transformers, loss-load factor and economic conductor sizing.
 *
 * References
 *  - Grainger & Stevenson, Power System Analysis (1994), ch. 4–6 (inductance from GMR/GMD, ABCD constants of the
 *    nominal-π and long-line models, surge-impedance loading).
 *  - Westinghouse Electrical Transmission and Distribution Reference Book (1964), ch. 4 (sequence impedances of
 *    cables with sheath and earth return, simplified Carson earth-return terms).
 *  - Carson (1926), Bell Syst. Tech. J. 5 (earth-return impedance; first-term approximation used here).
 *  - Kelvin's law for the economic conductor section (e.g. IEC 60287-3-2 for the full method).
 *
 * Phasors are rms, line-to-neutral voltages and line currents unless stated; complex values are { re, im }.
 * SI units: Ω/m, H/m, F/m, S/m, V, A, W, var, VA, m.
 */

import {
  cx, cadd, csub, cmul, cdiv, cscale, cconj, cabs, carg, csqrt, csinh, ccosh, abcdMul, toCx,
} from './complex.js';

const MU0 = 4e-7 * Math.PI;
const EPS0 = 8.8541878128e-12;
const TWO_PI = 2 * Math.PI;

function pos(v, name) {
  if (!(v > 0) || !Number.isFinite(v)) throw new RangeError(`${name}: must be a positive finite number`);
}
function nonNeg(v, name) {
  if (!(v >= 0) || !Number.isFinite(v)) throw new RangeError(`${name}: must be a non-negative finite number`);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Line constants                                                                                               */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Geometric mean distance between the phases of a three-phase circuit.
 * Trefoil (and the cores of a three-core cable): GMD = s; flat with equal spacing s: GMD = ∛2·s.
 * @param {object} o
 * @param {'trefoil'|'flat'|'threeCore'} o.formation
 * @param {number} o.s axial spacing between adjacent phases, m
 * @returns {number} GMD, m
 */
export function geometricMeanDistance({ formation, s }) {
  pos(s, 's');
  if (formation === 'trefoil' || formation === 'threeCore') return s;
  if (formation === 'flat') return Math.cbrt(2) * s;
  throw new RangeError('formation: must be trefoil, flat or threeCore');
}

/**
 * Positive-sequence inductance per phase, L = μ0/(2π)·ln(GMD/GMR), with GMR = gmrFactor·(dc/2).
 * gmrFactor = e^(−1/4) = 0.7788 for a solid round conductor with uniform current density; stranded conductors
 * have slightly lower factors (0.726 for 7 wires, 0.758 for 19, 0.768 for 37, 0.772 for 61, 0.774 for 91).
 * Reference: Grainger & Stevenson §4.
 * @param {object} o
 * @param {number} o.dc conductor diameter, m
 * @param {number} o.gmd geometric mean distance between phases, m
 * @param {number} [o.gmrFactor=0.7788]
 * @returns {number} inductance, H/m
 */
export function inductance({ dc, gmd, gmrFactor = Math.exp(-0.25) }) {
  pos(dc, 'dc'); pos(gmd, 'gmd'); pos(gmrFactor, 'gmrFactor');
  const gmr = gmrFactor * dc / 2;
  if (gmd <= gmr) throw new RangeError('gmd: must exceed the conductor GMR');
  return (MU0 / TWO_PI) * Math.log(gmd / gmr);
}

/**
 * Capacitance per phase of a screened core, C = 2π·ε0·εr/ln(Di/dc) (coaxial field; exact).
 * @param {object} o
 * @param {number} o.epsR relative permittivity
 * @param {number} o.Di diameter over the insulation, m
 * @param {number} o.dc diameter over the conductor screen, m
 * @returns {number} capacitance, F/m
 */
export function coaxialCapacitance({ epsR, Di, dc }) {
  pos(epsR, 'epsR'); pos(Di, 'Di'); pos(dc, 'dc');
  if (Di <= dc) throw new RangeError('Di: must exceed dc');
  return (TWO_PI * EPS0 * epsR) / Math.log(Di / dc);
}

/**
 * Charging current and reactive power of a cable (lumped, voltage assumed uniform):
 * Ic = ω·C·U0 per metre, Qc = 3·ω·C·U0² = ω·C·U² per metre of three-phase circuit.
 * @param {object} o
 * @param {number} o.C capacitance per phase, F/m
 * @param {number} o.U line-to-line voltage, V
 * @param {number} o.f frequency, Hz
 * @param {number} [o.length=1] m
 * @returns {{Ic:number, Qc:number, IcPerMetre:number, QcPerMetre:number}} A, var, A/m, var/m
 */
export function chargingCurrent({ C, U, f, length = 1 }) {
  nonNeg(C, 'C'); pos(U, 'U'); pos(f, 'f'); nonNeg(length, 'length');
  const w = TWO_PI * f;
  const IcPerMetre = (w * C * U) / Math.sqrt(3);
  const QcPerMetre = w * C * U * U;
  return { Ic: IcPerMetre * length, Qc: QcPerMetre * length, IcPerMetre, QcPerMetre };
}

/**
 * Critical length: the length at which the charging current of an uncompensated cable fed from one end equals
 * its current rating, l = I/(ω·C·U0). With equal compensation at both ends the usable length doubles.
 * @param {object} o
 * @param {number} o.C capacitance per phase, F/m
 * @param {number} o.U line-to-line voltage, V
 * @param {number} o.f frequency, Hz
 * @param {number} o.I current rating, A
 * @returns {number} critical length, m
 */
export function criticalLength({ C, U, f, I }) {
  pos(C, 'C'); pos(U, 'U'); pos(f, 'f'); pos(I, 'I');
  return (I * Math.sqrt(3)) / (TWO_PI * f * C * U);
}

/**
 * Surge (characteristic) impedance and natural loading of a lossless line: Zc = sqrt(L/C), SIL = U²/Zc,
 * propagation velocity v = 1/sqrt(L·C) and wavelength v/f.
 * Reference: Grainger & Stevenson §6.
 * @param {object} o
 * @param {number} o.L inductance per phase, H/m
 * @param {number} o.C capacitance per phase, F/m
 * @param {number} [o.U] line-to-line voltage for the SIL, V
 * @param {number} [o.f] frequency for the wavelength, Hz
 * @returns {{Zc:number, SIL:number, velocity:number, wavelength:number}} Ω, W, m/s, m (NaN when U or f missing)
 */
export function surgeImpedance({ L, C, U, f }) {
  pos(L, 'L'); pos(C, 'C');
  const Zc = Math.sqrt(L / C);
  const velocity = 1 / Math.sqrt(L * C);
  return { Zc, SIL: U ? (U * U) / Zc : NaN, velocity, wavelength: f ? velocity / f : NaN };
}

/**
 * Carson earth-return terms (first-term approximation): equivalent return depth De = 658.87·sqrt(ρ/f) m and
 * earth-return resistance Re = ω·μ0/8 = π²·f·10⁻⁷ Ω/m.
 * For a submarine cable the return path is the sea and seabed; use the seawater resistivity (about 0.2–0.3 Ω·m)
 * or an effective value for sea plus sediment. Valid at power frequency for conductor separations ≪ De.
 * Reference: Carson (1926); Westinghouse T&D Reference Book ch. 3.
 * @param {object} o
 * @param {number} o.f frequency, Hz
 * @param {number} o.rhoEarth resistivity of the return medium, Ω·m
 * @returns {{De:number, Re:number}} m, Ω/m
 */
export function earthReturn({ f, rhoEarth }) {
  pos(f, 'f'); pos(rhoEarth, 'rhoEarth');
  return { De: 658.87 * Math.sqrt(rhoEarth / f), Re: Math.PI * Math.PI * f * 1e-7 };
}

/**
 * Positive- and zero-sequence series impedance per unit length of a three-phase cable circuit made of three
 * screened cores (three single-core cables or a three-core cable with individually sheathed cores), with the
 * simplified Carson earth (sea) return.
 *
 *   zc = Rc + jωμ0/(2π)·ln(GMD/GMRc)                       conductor, positive sequence, sheaths open
 *   xm = ωμ0/(2π)·ln(GMD/rs)                               conductor–sheath mutual reactance
 *   Z1 = zc + xm²/(Rs + j·xm)                              sheaths bonded at both ends
 *   zc0 = Rc + 3Re + j·3ωμ0/(2π)·ln(De/GMR3c),  GMR3c = ∛(GMRc·GMD²)
 *   zs0 = Rs + 3Re + j·3ωμ0/(2π)·ln(De/GMR3s),  GMR3s = ∛(rs·GMD²)
 *   zm0 = 3Re + j·3ωμ0/(2π)·ln(De/GMR3s)
 *   Z0 = zc0 − zm0²/zs0  (return in sheaths and earth/sea in parallel)
 *   Z0 = zc0 + zs0 − 2·zm0  (return in the sheaths only)
 *   Z0 = zc0  (return in the earth/sea only, sheaths open)
 * Armour is not modelled separately: include it by using the parallel resistance of sheath and armour as Rs.
 * Reference: Westinghouse T&D Reference Book ch. 4 (cable sequence impedances).
 * @param {object} o
 * @param {number} o.Rc conductor AC resistance, Ω/m
 * @param {number} o.Rs sheath (screen) resistance, Ω/m
 * @param {number} o.dc conductor diameter, m
 * @param {number} o.ds mean sheath diameter, m
 * @param {number} o.gmd geometric mean distance between phases, m
 * @param {number} o.f frequency, Hz
 * @param {number} [o.rhoEarth=1] resistivity of the return medium, Ω·m (seawater ≈ 0.25, wet soil 10–100)
 * @param {number} [o.gmrFactor=0.7788] conductor GMR / radius
 * @param {'both'|'sheath'|'earth'} [o.zeroSequenceReturn='both']
 * @param {boolean} [o.sheathsBonded=true] sheaths bonded at both ends (circulating currents) for Z1
 * @returns {{Z1:{re,im}, Z0:{re,im}, zc0:{re,im}, zs0:{re,im}, zm0:{re,im}, De:number}} Ω/m; De in m
 */
export function sequenceImpedance({ Rc, Rs, dc, ds, gmd, f, rhoEarth = 1, gmrFactor = Math.exp(-0.25), zeroSequenceReturn = 'both', sheathsBonded = true }) {
  pos(Rc, 'Rc'); pos(Rs, 'Rs'); pos(dc, 'dc'); pos(ds, 'ds'); pos(gmd, 'gmd'); pos(f, 'f');
  const w = TWO_PI * f;
  const k = (w * MU0) / TWO_PI;
  const gmrC = gmrFactor * dc / 2;
  const rs = ds / 2;
  const { De, Re } = earthReturn({ f, rhoEarth });
  const zc = cx(Rc, k * Math.log(gmd / gmrC));
  const xm = k * Math.log(gmd / rs);
  const Z1 = sheathsBonded ? cadd(zc, cdiv(cx(xm * xm, 0), cx(Rs, xm))) : zc;
  const gmr3c = Math.cbrt(gmrC * gmd * gmd);
  const gmr3s = Math.cbrt(rs * gmd * gmd);
  const zc0 = cx(Rc + 3 * Re, 3 * k * Math.log(De / gmr3c));
  const zs0 = cx(Rs + 3 * Re, 3 * k * Math.log(De / gmr3s));
  const zm0 = cx(3 * Re, 3 * k * Math.log(De / gmr3s));
  let Z0;
  if (zeroSequenceReturn === 'both') Z0 = csub(zc0, cdiv(cmul(zm0, zm0), zs0));
  else if (zeroSequenceReturn === 'sheath') Z0 = csub(cadd(zc0, zs0), cscale(zm0, 2));
  else if (zeroSequenceReturn === 'earth') Z0 = zc0;
  else throw new RangeError('zeroSequenceReturn: must be both, sheath or earth');
  return { Z1, Z0, zc0, zs0, zm0, De };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Two-port line models                                                                                         */
/* ------------------------------------------------------------------------------------------------------------ */

function lineZY({ r, l, c, g = 0, f }) {
  nonNeg(r, 'r'); nonNeg(l, 'l'); nonNeg(c, 'c'); nonNeg(g, 'g'); pos(f, 'f');
  const w = TWO_PI * f;
  return { z: cx(r, w * l), y: cx(g, w * c) };
}

/**
 * Propagation constant and characteristic impedance of a line: γ = sqrt(z·y), Zc = sqrt(z/y),
 * z = r + jωl, y = g + jωc.
 * @param {object} o
 * @param {number} o.r Ω/m @param {number} o.l H/m @param {number} o.c F/m @param {number} [o.g=0] S/m
 * @param {number} o.f frequency, Hz
 * @returns {{gamma:{re,im}, Zc:{re,im}, z:{re,im}, y:{re,im}}} 1/m, Ω, Ω/m, S/m
 */
export function propagation({ r, l, c, g = 0, f }) {
  const { z, y } = lineZY({ r, l, c, g, f });
  if (cabs(y) === 0 || cabs(z) === 0) throw new RangeError('l, c: line needs non-zero series impedance and shunt admittance');
  return { gamma: csqrt(cmul(z, y)), Zc: csqrt(cdiv(z, y)), z, y };
}

/**
 * ABCD constants of the exact distributed-parameter line: A = D = cosh(γl), B = Zc·sinh(γl), C = sinh(γl)/Zc.
 * Sending-end quantities follow from Vs = A·Vr + B·Ir, Is = C·Vr + D·Ir.
 * Reference: Grainger & Stevenson §6.
 * @param {object} o line constants r, l, c, g, f (see propagation) and
 * @param {number} o.length m
 * @returns {{A:{re,im}, B:{re,im}, C:{re,im}, D:{re,im}}}
 */
export function abcdDistributed({ r, l, c, g = 0, f, length }) {
  nonNeg(length, 'length');
  const { gamma, Zc } = propagation({ r, l, c, g, f });
  const gl = cscale(gamma, length);
  const ch = ccosh(gl); const sh = csinh(gl);
  return { A: ch, B: cmul(Zc, sh), C: cdiv(sh, Zc), D: ch };
}

/**
 * ABCD constants of the nominal-π model: A = D = 1 + Z·Y/2, B = Z, C = Y·(1 + Z·Y/4) with Z = z·l, Y = y·l.
 * Adequate for cables up to a few tens of km at 50/60 Hz; use the distributed model beyond.
 * @param {object} o line constants r, l, c, g, f and
 * @param {number} o.length m
 * @returns {{A:{re,im}, B:{re,im}, C:{re,im}, D:{re,im}}}
 */
export function abcdNominalPi({ r, l, c, g = 0, f, length }) {
  nonNeg(length, 'length');
  const { z, y } = lineZY({ r, l, c, g, f });
  const Z = cscale(z, length); const Y = cscale(y, length);
  const zy = cmul(Z, Y);
  const A = cadd(cx(1, 0), cscale(zy, 0.5));
  return { A, B: Z, C: cmul(Y, cadd(cx(1, 0), cscale(zy, 0.25))), D: A };
}

/**
 * ABCD constants of a shunt admittance (reactor, capacitor bank): [[1, 0], [Y, 1]].
 * @param {{re:number, im:number}|number} Y shunt admittance per phase, S (inductive reactor: negative imaginary)
 * @returns {{A, B, C, D}}
 */
export function abcdShunt(Y) {
  return { A: cx(1, 0), B: cx(0, 0), C: toCx(Y), D: cx(1, 0) };
}

/**
 * ABCD constants of a series impedance: [[1, Z], [0, 1]].
 * @param {{re:number, im:number}|number} Z series impedance per phase, Ω
 * @returns {{A, B, C, D}}
 */
export function abcdSeries(Z) {
  return { A: cx(1, 0), B: toCx(Z), C: cx(0, 0), D: cx(1, 0) };
}

/**
 * Cascade of two-ports from the sending end to the receiving end.
 * @param {{A,B,C,D}[]} parts two-ports in order sending → receiving
 * @returns {{A, B, C, D}}
 */
export function abcdCascade(parts) {
  if (!Array.isArray(parts) || parts.length === 0) throw new RangeError('parts: must be a non-empty array');
  return parts.reduce((m, p) => abcdMul(m, p));
}

/**
 * Operating point of a three-phase two-port given the conditions at one end.
 * Receiving end known: Vs = A·Vr + B·Ir, Is = C·Vr + D·Ir. Sending end known: Vr = D·Vs − B·Is, Ir = −C·Vs + A·Is
 * (reciprocal network, AD − BC = 1).
 * P and Q are three-phase totals flowing from the sending towards the receiving end at the known end; Q > 0 is
 * lagging (inductive) power flowing in that direction.
 * @param {object} o
 * @param {{A,B,C,D}} o.abcd two-port constants
 * @param {'receiving'|'sending'} [o.known='receiving'] end at which U, P, Q are specified
 * @param {number} o.U line-to-line voltage at the known end, V
 * @param {number} o.P active power at the known end, W
 * @param {number} [o.Q=0] reactive power at the known end, var
 * @returns {{Vs:{re,im}, Is:{re,im}, Vr:{re,im}, Ir:{re,im}, Us:number, Ur:number, Ps:number, Qs:number, Pr:number,
 *            Qr:number, losses:number, efficiency:number, regulation:number, voltageDrop:number, angle:number}}
 *          phasors (line-to-neutral V, A), line-to-line magnitudes (V), three-phase powers (W, var), losses (W),
 *          efficiency Pr/Ps, regulation (|Vr,no-load| − |Vr|)/|Vr| with Vr,no-load = Vs/A, voltage drop
 *          (Us − Ur)/Ur and transmission angle δ = arg(Vs) − arg(Vr) (rad)
 */
export function solveTwoPort({ abcd, known = 'receiving', U, P, Q = 0 }) {
  pos(U, 'U');
  if (!Number.isFinite(P) || !Number.isFinite(Q)) throw new RangeError('P, Q: must be finite');
  const V = cx(U / Math.sqrt(3), 0);
  // S = 3·V·I*  ⇒  I = conj(S/(3V))
  const I = cconj(cdiv(cx(P, Q), cscale(V, 3)));
  let Vs; let Is; let Vr; let Ir;
  if (known === 'receiving') {
    Vr = V; Ir = I;
    Vs = cadd(cmul(abcd.A, Vr), cmul(abcd.B, Ir));
    Is = cadd(cmul(abcd.C, Vr), cmul(abcd.D, Ir));
  } else if (known === 'sending') {
    Vs = V; Is = I;
    Vr = csub(cmul(abcd.D, Vs), cmul(abcd.B, Is));
    Ir = csub(cmul(abcd.A, Is), cmul(abcd.C, Vs));
  } else throw new RangeError('known: must be receiving or sending');
  const Ss = cscale(cmul(Vs, cconj(Is)), 3);
  const Sr = cscale(cmul(Vr, cconj(Ir)), 3);
  const Us = cabs(Vs) * Math.sqrt(3); const Ur = cabs(Vr) * Math.sqrt(3);
  const noLoad = cabs(cdiv(Vs, abcd.A));
  return { Vs, Is, Vr, Ir, Us, Ur, Ps: Ss.re, Qs: Ss.im, Pr: Sr.re, Qr: Sr.im, losses: Ss.re - Sr.re,
    efficiency: Ss.re !== 0 ? Sr.re / Ss.re : NaN, regulation: (noLoad - cabs(Vr)) / cabs(Vr),
    voltageDrop: (Us - Ur) / Ur, angle: carg(Vs) - carg(Vr) };
}

/**
 * Voltage and current magnitude along a uniform line from its receiving-end phasors:
 * V(x) = cosh(γx)·Vr + Zc·sinh(γx)·Ir, I(x) = sinh(γx)/Zc·Vr + cosh(γx)·Ir, x measured from the receiving end.
 * In a cable the current is normally largest at one of the ends (charging current adds to the load current).
 * @param {object} o line constants r, l, c, g, f, length and
 * @param {{re:number, im:number}} o.Vr receiving-end voltage (line-to-neutral), V
 * @param {{re:number, im:number}} o.Ir receiving-end current, A
 * @param {number} [o.n=51] number of points
 * @returns {{x:Float64Array, U:Float64Array, I:Float64Array, Imax:number, xImax:number, Umax:number}}
 *          x from the receiving end (m), line-to-line voltage (V), current (A), maximum current and its position
 */
export function lineProfile({ r, l, c, g = 0, f, length, Vr, Ir, n = 51 }) {
  pos(length, 'length');
  if (!(n >= 2)) throw new RangeError('n: at least 2 points');
  const { gamma, Zc } = propagation({ r, l, c, g, f });
  const x = new Float64Array(n); const U = new Float64Array(n); const I = new Float64Array(n);
  let Imax = 0; let xImax = 0; let Umax = 0;
  for (let i = 0; i < n; i++) {
    x[i] = (length * i) / (n - 1);
    const gx = cscale(gamma, x[i]);
    const ch = ccosh(gx); const sh = csinh(gx);
    const V = cadd(cmul(ch, Vr), cmul(cmul(Zc, sh), Ir));
    const Ii = cadd(cmul(cdiv(sh, Zc), Vr), cmul(ch, Ir));
    U[i] = cabs(V) * Math.sqrt(3); I[i] = cabs(Ii);
    if (I[i] > Imax) { Imax = I[i]; xImax = x[i]; }
    if (U[i] > Umax) Umax = U[i];
  }
  return { x, U, I, Imax, xImax, Umax };
}

/**
 * Open-circuit (no-load) voltage rise of a line energised from one end (Ferranti effect): Vr/Vs = 1/cosh(γl).
 * @param {object} o line constants r, l, c, g, f and length (m)
 * @returns {{ratio:number, rise:number}} |Vr/Vs| and the rise |Vr/Vs| − 1
 */
export function ferrantiRise({ r, l, c, g = 0, f, length }) {
  const { gamma } = propagation({ r, l, c, g, f });
  const ratio = 1 / cabs(ccosh(cscale(gamma, length)));
  return { ratio, rise: ratio - 1 };
}

/**
 * Compensated cable link: sending bus → shunt reactor → cable (→ mid-point reactor → cable) → shunt reactor →
 * receiving bus. Reactor sizes are fractions of the total cable charging power ω·C·U²·l at nominal voltage.
 * The sending end (e.g. offshore substation) is the reference: its voltage U, active power P and reactive power Q
 * (as seen by the link, downstream of nothing) are given; the cable currents at every section end and the
 * receiving-end quantities follow from the exact distributed model.
 * @param {object} o line constants r, l, c, g, f, length and
 * @param {number} o.U sending-end line-to-line voltage, V
 * @param {number} o.P sending-end active power, W
 * @param {number} [o.Q=0] sending-end reactive power into the link, var
 * @param {{sending?:number, mid?:number, receiving?:number}} [o.compensation] reactor fractions (default none)
 * @param {number} [o.nProfile=21] points per cable section for the current profile
 * @returns {{Pr:number, Qr:number, Ur:number, losses:number, efficiency:number, ImaxCable:number, Is:number, Ir:number,
 *            reactors:{sending:number, mid:number, receiving:number}, Qcharging:number}}
 *          receiving power (W, var), receiving voltage (V), losses (W), largest cable current (A), cable current at
 *          the sending and receiving ends (A), reactor ratings at nominal voltage (var) and total charging power (var)
 */
export function compensatedLink({ r, l, c, g = 0, f, length, U, P, Q = 0, compensation = {}, nProfile = 21 }) {
  pos(length, 'length'); pos(U, 'U');
  const kS = compensation.sending ?? 0; const kM = compensation.mid ?? 0; const kR = compensation.receiving ?? 0;
  nonNeg(kS, 'compensation.sending'); nonNeg(kM, 'compensation.mid'); nonNeg(kR, 'compensation.receiving');
  const w = TWO_PI * f;
  const Bc = w * c * length; // total cable susceptance per phase, S
  const shunt = (k) => abcdShunt(cx(0, -k * Bc));
  const V = cx(U / Math.sqrt(3), 0);
  let Vx = V;
  let Ix = cconj(cdiv(cx(P, Q), cscale(V, 3)));
  const step = (m) => {
    // Advance from the sending side of a two-port to its receiving side.
    const Vn = csub(cmul(m.D, Vx), cmul(m.B, Ix));
    const In = csub(cmul(m.A, Ix), cmul(m.C, Vx));
    Vx = Vn; Ix = In;
  };
  let Imax = 0;
  const cableSection = (len) => {
    // Largest current inside the section: sample the exact profile from its receiving end.
    step(abcdDistributed({ r, l, c, g, f, length: len }));
    const prof = lineProfile({ r, l, c, g, f, length: len, Vr: Vx, Ir: Ix, n: nProfile });
    if (prof.Imax > Imax) Imax = prof.Imax;
  };
  if (kS > 0) step(shunt(kS));
  const IsCable = cabs(Ix);
  if (kM > 0) { cableSection(length / 2); step(shunt(kM)); cableSection(length / 2); } else cableSection(length);
  const IrCable = cabs(Ix);
  if (kR > 0) step(shunt(kR));
  const Sr = cscale(cmul(Vx, cconj(Ix)), 3);
  const Qch = w * c * U * U * length;
  return { Pr: Sr.re, Qr: Sr.im, Ur: cabs(Vx) * Math.sqrt(3), losses: P - Sr.re, efficiency: P !== 0 ? Sr.re / P : NaN,
    ImaxCable: Imax, Is: IsCable, Ir: IrCable,
    reactors: { sending: kS * Qch, mid: kM * Qch, receiving: kR * Qch }, Qcharging: Qch };
}

/**
 * Active-power transfer capability of an AC cable link as a function of length: for each length the largest
 * sending-end active power for which the cable current nowhere exceeds the rating (bisection on P with the
 * exact distributed model, sending end at nominal voltage and the stated reactive power).
 * @param {object} o line constants r, l, c, g, f and
 * @param {number[]} o.lengths m
 * @param {number} o.U line-to-line voltage, V
 * @param {number} o.Irated cable current rating, A
 * @param {number} [o.Q=0] sending-end reactive power into the link, var
 * @param {{sending?:number, mid?:number, receiving?:number}} [o.compensation]
 * @returns {{lengths:number[], Psend:number[], Preceive:number[], efficiency:number[]}} W; zeros where the
 *          charging current alone exceeds the rating
 */
export function transferCapability({ r, l, c, g = 0, f, lengths, U, Irated, Q = 0, compensation = {} }) {
  pos(Irated, 'Irated'); pos(U, 'U');
  const Psend = []; const Preceive = []; const efficiency = [];
  const Smax = Math.sqrt(3) * U * Irated;
  for (const length of lengths) {
    const run = (P) => compensatedLink({ r, l, c, g, f, length, U, P, Q, compensation });
    if (run(0).ImaxCable >= Irated) { Psend.push(0); Preceive.push(0); efficiency.push(0); continue; }
    let lo = 0; let hi = 1.5 * Smax;
    for (let i = 0; i < 50; i++) { const mid = 0.5 * (lo + hi); if (run(mid).ImaxCable > Irated) hi = mid; else lo = mid; }
    const res = run(lo);
    Psend.push(lo); Preceive.push(Math.max(0, res.Pr)); efficiency.push(lo > 0 ? Math.max(0, res.Pr) / lo : 0);
  }
  return { lengths: Array.from(lengths), Psend, Preceive, efficiency };
}

/**
 * Split of a given total shunt compensation between the two ends of a cable that maximises the transferable
 * active power (scan of the sending-end share; for a symmetrical link the optimum is close to one half).
 * @param {object} o as transferCapability with a single `length` (m) and
 * @param {number} [o.total=1] total compensation as a fraction of the cable charging power
 * @param {number} [o.steps=20] resolution of the scan
 * @returns {{sendingShare:number, Preceive:number, table:{share:number, Preceive:number}[]}}
 */
export function optimalCompensationSplit({ r, l, c, g = 0, f, length, U, Irated, Q = 0, total = 1, steps = 20 }) {
  nonNeg(total, 'total');
  if (!(steps >= 2)) throw new RangeError('steps: at least 2');
  const table = []; let best = { share: 0, Preceive: -1 };
  for (let i = 0; i <= steps; i++) {
    const share = i / steps;
    const cap = transferCapability({ r, l, c, g, f, lengths: [length], U, Irated, Q,
      compensation: { sending: share * total, receiving: (1 - share) * total } });
    const row = { share, Preceive: cap.Preceive[0] };
    table.push(row);
    if (row.Preceive > best.Preceive) best = row;
  }
  return { sendingShare: best.share, Preceive: best.Preceive, table };
}

/**
 * Resonance estimates of a cable connected to an inductive grid.
 * Parallel (first) resonance of the lumped cable capacitance with the source inductance,
 * f = 1/(2π·sqrt(Lgrid·C·l)); quarter-wave resonance of the open-ended cable itself, f = 1/(4·l·sqrt(L·C));
 * and the first resonance of the distributed lossless cable behind the grid inductance, the lowest root of
 * ω·Lgrid·tan(ω·l·sqrt(L·C)) = Zc (input admittance zero), which tends to the two estimates in the limits.
 * @param {object} o
 * @param {number} o.L cable inductance, H/m @param {number} o.C cable capacitance, F/m @param {number} o.length m
 * @param {number} o.Lgrid source inductance per phase, H (= U²/(2π·f·S″k))
 * @returns {{fLumped:number, fQuarterWave:number, fDistributed:number}} Hz
 */
export function resonanceFrequency({ L, C, length, Lgrid }) {
  pos(L, 'L'); pos(C, 'C'); pos(length, 'length'); pos(Lgrid, 'Lgrid');
  const fLumped = 1 / (TWO_PI * Math.sqrt(Lgrid * C * length));
  const tau = length * Math.sqrt(L * C);
  const fQuarterWave = 1 / (4 * tau);
  const Zc = Math.sqrt(L / C);
  // Root of w·Lgrid·tan(w·tau) − Zc on (0, π/(2·tau)): monotonic, bracketed.
  let lo = 0; let hi = (Math.PI / (2 * tau)) * (1 - 1e-12);
  for (let i = 0; i < 200; i++) {
    const mid = 0.5 * (lo + hi);
    if (mid * Lgrid * Math.tan(mid * tau) < Zc) lo = mid; else hi = mid;
  }
  return { fLumped, fQuarterWave, fDistributed: (0.5 * (lo + hi)) / TWO_PI };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* HVDC                                                                                                         */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Steady-state performance of an HVDC cable link.
 * Loop resistance: monopole with sea/earth return Rloop = r·l + Relectrodes; monopole with metallic return
 * Rloop = (r + rReturn)·l; symmetric monopole and bipole (balanced) Rloop = 2·r·l with pole-to-pole voltage 2·U.
 * P = Upp·I at the sending converter DC terminals; cable loss I²·Rloop; each converter station loses a fraction
 * of the power passing through it (typical 0.7–1.0 % per station for VSC, 0.7–0.8 % for LCC).
 * @param {object} o
 * @param {'monopoleSeaReturn'|'monopoleMetallicReturn'|'symmetricMonopole'|'bipole'} o.configuration
 * @param {number} o.U pole voltage to earth, V
 * @param {number} o.I pole current, A
 * @param {number} o.r DC resistance of one pole conductor at operating temperature, Ω/m
 * @param {number} o.length route length, m
 * @param {number} [o.rReturn=r] resistance of the metallic return conductor, Ω/m
 * @param {number} [o.electrodeResistance=0] total electrode resistance for sea return, Ω
 * @param {number} [o.converterLossFraction=0] loss per converter station as a fraction of throughput
 * @returns {{Pdc:number, Rloop:number, cableLoss:number, voltageDrop:number, voltageDropFraction:number,
 *            converterLoss:number, Pin:number, Pout:number, efficiency:number}} W, Ω, W, V, –, W, W, W, –
 *          (Pin = AC power into the sending converter, Pout = AC power out of the receiving converter;
 *          voltageDrop is pole-to-pole or pole-to-return across the loop)
 */
export function hvdcLink({ configuration, U, I, r, length, rReturn, electrodeResistance = 0, converterLossFraction = 0 }) {
  pos(U, 'U'); nonNeg(I, 'I'); nonNeg(r, 'r'); nonNeg(length, 'length');
  if (!(converterLossFraction >= 0 && converterLossFraction < 1)) throw new RangeError('converterLossFraction: must be in [0, 1)');
  let Rloop; let Upp;
  switch (configuration) {
    case 'monopoleSeaReturn': Rloop = r * length + electrodeResistance; Upp = U; break;
    case 'monopoleMetallicReturn': Rloop = (r + (rReturn ?? r)) * length; Upp = U; break;
    case 'symmetricMonopole': case 'bipole': Rloop = 2 * r * length; Upp = 2 * U; break;
    default: throw new RangeError('configuration: must be monopoleSeaReturn, monopoleMetallicReturn, symmetricMonopole or bipole');
  }
  const Pdc = Upp * I;
  const cableLoss = I * I * Rloop;
  const voltageDrop = I * Rloop;
  const Pin = Pdc / (1 - converterLossFraction);
  const PdcReceive = Pdc - cableLoss;
  const Pout = PdcReceive * (1 - converterLossFraction);
  return { Pdc, Rloop, cableLoss, voltageDrop, voltageDropFraction: voltageDrop / Upp,
    converterLoss: (Pin - Pdc) + (PdcReceive - Pout), Pin, Pout, efficiency: Pin > 0 ? Pout / Pin : NaN };
}

/**
 * Delivered power and loss fraction of an HVAC link and an HVDC link over a range of lengths, and the length
 * beyond which the DC link delivers more power (or, with cost functions supplied by the caller, costs less per
 * delivered watt). No cost data are built in.
 * @param {object} o
 * @param {number[]} o.lengths route lengths, m (ascending)
 * @param {object} o.ac arguments of transferCapability without `lengths` (r, l, c, f, U, Irated, compensation ...)
 * @param {object} o.dc arguments of hvdcLink without `length` (configuration, U, I, r, converterLossFraction ...)
 * @param {(length:number)=>number} [o.costAc] total cost of the AC link for a length (any currency)
 * @param {(length:number)=>number} [o.costDc] total cost of the DC link
 * @returns {{lengths:number[], Pac:number[], Pdc:number[], lossFractionAc:number[], lossFractionDc:number[],
 *            breakEvenCapacity:?number, breakEvenCost:?number}} W, –, and break-even lengths in m (null if the
 *          curves do not cross in the range). breakEvenCost compares cost per delivered watt.
 */
export function acDcBreakEven({ lengths, ac, dc, costAc, costDc }) {
  if (!Array.isArray(lengths) || lengths.length < 2) throw new RangeError('lengths: at least two values');
  const cap = transferCapability({ ...ac, lengths });
  const Pac = cap.Preceive; const Pdc = []; const lossFractionAc = []; const lossFractionDc = [];
  for (let i = 0; i < lengths.length; i++) {
    const h = hvdcLink({ ...dc, length: lengths[i] });
    Pdc.push(h.Pout);
    lossFractionDc.push(1 - h.efficiency);
    lossFractionAc.push(cap.Psend[i] > 0 ? 1 - cap.efficiency[i] : 1);
  }
  const cross = (fa, fb) => {
    // First length at which fb − fa changes sign from negative to non-negative (linear interpolation).
    for (let i = 1; i < lengths.length; i++) {
      const d0 = fb[i - 1] - fa[i - 1]; const d1 = fb[i] - fa[i];
      if (d0 < 0 && d1 >= 0) return lengths[i - 1] + ((lengths[i] - lengths[i - 1]) * -d0) / (d1 - d0);
    }
    return null;
  };
  const breakEvenCapacity = cross(Pac, Pdc);
  let breakEvenCost = null;
  if (costAc && costDc) {
    const perWattAc = lengths.map((L, i) => (Pac[i] > 0 ? costAc(L) / Pac[i] : Infinity));
    const perWattDc = lengths.map((L, i) => (Pdc[i] > 0 ? costDc(L) / Pdc[i] : Infinity));
    // DC becomes cheaper where (AC − DC) cost per watt turns non-negative.
    const a = perWattDc.map((v) => (Number.isFinite(v) ? v : 1e300));
    const b = perWattAc.map((v) => (Number.isFinite(v) ? v : 1e300));
    breakEvenCost = cross(a, b);
  }
  return { lengths: Array.from(lengths), Pac, Pdc, lossFractionAc, lossFractionDc, breakEvenCapacity, breakEvenCost };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Transformers, per-unit, three-phase relations                                                                */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Per-unit base quantities of a three-phase system: Zbase = U²/S, Ibase = S/(√3·U).
 * @param {object} o
 * @param {number} o.S base apparent power (three-phase), VA
 * @param {number} o.U base line-to-line voltage, V
 * @returns {{Zbase:number, Ibase:number, Ybase:number}} Ω, A, S
 */
export function perUnitBase({ S, U }) {
  pos(S, 'S'); pos(U, 'U');
  const Zbase = (U * U) / S;
  return { Zbase, Ibase: S / (Math.sqrt(3) * U), Ybase: 1 / Zbase };
}

/**
 * Convert a per-unit impedance between bases: Znew = Zold·(Uold/Unew)²·(Snew/Sold).
 * @param {object} o
 * @param {number} o.z per-unit value on the old base
 * @param {number} o.Sold VA @param {number} o.Uold V @param {number} o.Snew VA @param {number} o.Unew V
 * @returns {number} per-unit value on the new base
 */
export function perUnitChangeBase({ z, Sold, Uold, Snew, Unew }) {
  pos(Sold, 'Sold'); pos(Uold, 'Uold'); pos(Snew, 'Snew'); pos(Unew, 'Unew');
  return z * (Uold / Unew) ** 2 * (Snew / Sold);
}

/**
 * Balanced three-phase power relations: S = √3·U·I, P = S·cosφ, Q = S·sinφ (Q > 0 lagging).
 * Give any two of (P, S, I) together with U and the power factor.
 * @param {object} o
 * @param {number} o.U line-to-line voltage, V
 * @param {number} [o.I] line current, A @param {number} [o.P] active power, W @param {number} [o.S] apparent power, VA
 * @param {number} [o.powerFactor=1] cosφ (0–1)
 * @param {boolean} [o.leading=false] leading (capacitive) power factor
 * @returns {{S:number, P:number, Q:number, I:number, phi:number}} VA, W, var, A, rad
 */
export function threePhasePower({ U, I, P, S, powerFactor = 1, leading = false }) {
  pos(U, 'U');
  if (!(powerFactor > 0 && powerFactor <= 1)) throw new RangeError('powerFactor: must be in (0, 1]');
  let s;
  if (S !== undefined) s = S;
  else if (I !== undefined) s = Math.sqrt(3) * U * I;
  else if (P !== undefined) s = P / powerFactor;
  else throw new RangeError('I, P or S: one of them is required');
  const phi = Math.acos(powerFactor) * (leading ? -1 : 1);
  return { S: s, P: s * powerFactor, Q: s * Math.sin(phi), I: s / (Math.sqrt(3) * U), phi };
}

/**
 * Two-winding transformer from its test data: series impedance, regulation and efficiency at a load.
 * Per-unit series resistance er = Pk/Sr (load loss at rated current), impedance ez = uk, reactance
 * ex = sqrt(uk² − er²). Regulation is evaluated exactly from the phasor equation
 * |V1| = |1 + (er + j·ex)·x·(cosφ − j·sinφ)| (secondary at 1 pu), not from the two-term approximation.
 * Efficiency = x·Sr·cosφ/(x·Sr·cosφ + P0 + x²·Pk); it peaks at x = sqrt(P0/Pk).
 * Reference: IEC 60076-1 definitions; any machines text.
 * @param {object} o
 * @param {number} o.Sr rated power, VA
 * @param {number} o.uk short-circuit impedance voltage, per unit (e.g. 0.12)
 * @param {number} o.Pk load loss at rated current, W
 * @param {number} o.P0 no-load loss, W
 * @param {number} [o.loading=1] load as a fraction of rated apparent power x
 * @param {number} [o.powerFactor=1] load power factor @param {boolean} [o.leading=false]
 * @param {number} [o.U] rated line-to-line voltage of the winding to which Z is referred, V
 * @returns {{er:number, ex:number, Z:{re:number,im:number}|null, regulation:number, efficiency:number, losses:number,
 *            loadingAtMaxEfficiency:number, maxEfficiency:number}} per-unit er, ex; Z in Ω (null without U);
 *          regulation as a fraction; losses in W
 */
export function transformer({ Sr, uk, Pk, P0, loading = 1, powerFactor = 1, leading = false, U }) {
  pos(Sr, 'Sr'); pos(uk, 'uk'); nonNeg(Pk, 'Pk'); nonNeg(P0, 'P0'); nonNeg(loading, 'loading');
  if (!(powerFactor > 0 && powerFactor <= 1)) throw new RangeError('powerFactor: must be in (0, 1]');
  const er = Pk / Sr;
  if (er > uk) throw new RangeError('Pk: load loss implies a resistance larger than the impedance uk');
  const ex = Math.sqrt(uk * uk - er * er);
  const sinPhi = Math.sqrt(1 - powerFactor * powerFactor) * (leading ? -1 : 1);
  const i = cx(loading * powerFactor, -loading * sinPhi);
  const v1 = cadd(cx(1, 0), cmul(cx(er, ex), i));
  const Pout = loading * Sr * powerFactor;
  const losses = P0 + loading * loading * Pk;
  const xOpt = Pk > 0 ? Math.sqrt(P0 / Pk) : Infinity;
  const eff = (x) => (x * Sr * powerFactor) / (x * Sr * powerFactor + P0 + x * x * Pk);
  return { er, ex, Z: U ? cscale(cx(er, ex), (U * U) / Sr) : null, regulation: cabs(v1) - 1,
    efficiency: Pout > 0 ? Pout / (Pout + losses) : 0, losses,
    loadingAtMaxEfficiency: xOpt, maxEfficiency: Number.isFinite(xOpt) && xOpt > 0 ? eff(xOpt) : NaN };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Losses over time and economic sizing                                                                         */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Load factor and loss-load factor of a load profile or load-duration curve, and the annual energy loss.
 * LF = mean(I)/Imax, LLF = mean(I²)/Imax² (exact for current-dependent losses at constant resistance);
 * annual loss = lossAtPeak·LLF·hours. The empirical relation LLF ≈ k·LF + (1 − k)·LF² (k ≈ 0.1–0.3) is returned
 * for comparison only.
 * @param {object} o
 * @param {ArrayLike<number>} o.load load samples (current or power at constant voltage and power factor), any unit
 * @param {ArrayLike<number>} [o.weights] duration or probability of each sample (default equal)
 * @param {number} [o.lossAtPeak] loss at the peak load, W (e.g. 3·Imax²·R·l)
 * @param {number} [o.hours=8760] hours in the period
 * @param {number} [o.k=0.2] coefficient of the empirical relation
 * @returns {{loadFactor:number, lossLoadFactor:number, lossLoadFactorEmpirical:number, peak:number,
 *            energyLoss:number, equivalentLossHours:number}} energyLoss in Wh (NaN without lossAtPeak)
 */
export function lossLoadFactor({ load, weights, lossAtPeak, hours = 8760, k = 0.2 }) {
  const n = load.length;
  if (!(n >= 1)) throw new RangeError('load: must be a non-empty array');
  if (weights && weights.length !== n) throw new RangeError('weights: must have the same length as load');
  let peak = 0; let sw = 0; let s1 = 0; let s2 = 0;
  for (let i = 0; i < n; i++) {
    const v = Math.abs(load[i]); const w = weights ? weights[i] : 1;
    if (!(w >= 0)) throw new RangeError('weights: must be ≥ 0');
    if (v > peak) peak = v;
    sw += w; s1 += w * v; s2 += w * v * v;
  }
  if (!(peak > 0) || !(sw > 0)) throw new RangeError('load: peak and total weight must be positive');
  const loadFactor = s1 / sw / peak;
  const llf = s2 / sw / (peak * peak);
  return { loadFactor, lossLoadFactor: llf, lossLoadFactorEmpirical: k * loadFactor + (1 - k) * loadFactor * loadFactor,
    peak, energyLoss: lossAtPeak !== undefined ? lossAtPeak * llf * hours : NaN, equivalentLossHours: llf * hours };
}

/**
 * Capital recovery factor (annuity): crf = i·(1 + i)^N/((1 + i)^N − 1); 1/N for i = 0.
 * @param {object} o
 * @param {number} o.rate discount rate per year (fraction) @param {number} o.years lifetime, years
 * @returns {number} crf, 1/year
 */
export function capitalRecoveryFactor({ rate, years }) {
  nonNeg(rate, 'rate'); pos(years, 'years');
  if (rate === 0) return 1 / years;
  const a = Math.pow(1 + rate, years);
  return (rate * a) / (a - 1);
}

/**
 * Economic conductor cross-section (Kelvin's law). Annual cost per metre of circuit
 *   c(S) = crf·(a + b·S) + nPhases·Imax²·ρ/S·LLF·hours·energyPrice
 * is minimised at S* = Imax·sqrt(nPhases·ρ·LLF·hours·energyPrice/(crf·b)), where the annualised
 * section-dependent capital cost equals the annual cost of losses. With a list of candidate sections and their
 * installed costs the cheapest candidate is also returned. Sheath/armour/dielectric losses and load growth are
 * not included (see IEC 60287-3-2 for the complete method); all cost data come from the caller.
 * @param {object} o
 * @param {number} o.Imax peak current, A
 * @param {number} o.rho conductor AC resistivity at the mean operating temperature, Ω·m
 * @param {number} o.lossLoadFactor LLF of the load cycle (0–1)
 * @param {number} o.energyPrice value of lost energy, currency/Wh
 * @param {number} o.costPerArea section-dependent installed cost b, currency per metre of circuit per m² of
 *        conductor section
 * @param {number} [o.fixedCost=0] section-independent cost a, currency/m
 * @param {number} o.crf capital recovery factor, 1/year
 * @param {number} [o.nPhases=3] loaded conductors
 * @param {number} [o.hours=8760]
 * @param {{area:number, cost:number, R?:number}[]} [o.candidates] standard sections (m²) with installed cost
 *        (currency/m) and optionally their own AC resistance (Ω/m)
 * @returns {{areaOptimal:number, annualCostOptimal:number, currentDensityOptimal:number,
 *            best:?{area:number, annualCost:number, capital:number, losses:number}, table:object[]}}
 *          m², currency/(m·year), A/m²
 */
export function economicConductorSize({ Imax, rho, lossLoadFactor: llf, energyPrice, costPerArea, fixedCost = 0, crf,
  nPhases = 3, hours = 8760, candidates }) {
  pos(Imax, 'Imax'); pos(rho, 'rho'); pos(llf, 'lossLoadFactor'); pos(energyPrice, 'energyPrice'); pos(costPerArea, 'costPerArea');
  pos(crf, 'crf'); nonNeg(fixedCost, 'fixedCost');
  const kLoss = nPhases * Imax * Imax * llf * hours * energyPrice; // annual loss cost = kLoss·R
  const areaOptimal = Math.sqrt((kLoss * rho) / (crf * costPerArea));
  const annualCostOptimal = crf * (fixedCost + costPerArea * areaOptimal) + (kLoss * rho) / areaOptimal;
  let best = null; const table = [];
  if (candidates) {
    for (const cnd of candidates) {
      pos(cnd.area, 'candidates[].area');
      const capital = crf * cnd.cost;
      const losses = kLoss * (cnd.R ?? rho / cnd.area);
      const row = { area: cnd.area, annualCost: capital + losses, capital, losses };
      table.push(row);
      if (!best || row.annualCost < best.annualCost) best = row;
    }
  }
  return { areaOptimal, annualCostOptimal, currentDensityOptimal: Imax / areaOptimal, best, table };
}
