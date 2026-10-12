/**
 * Static strength of simple (non-overlapping, unstiffened, ungrouted) circular tubular joints to
 * API RP 2A-WSD 21st ed. Supplement 2/3 and 22nd ed., Section 4.3 (basic capacity, strength factor Qu,
 * chord load factor Qf, interaction and the minimum-capacity requirement of §4.2.3).
 *
 * Geometry parameters: β = d/D, γ = D/(2T), τ = t/T, θ = brace angle from the chord, g = gap.
 * SI units throughout (m, N, N·m, Pa, rad; arguments ending in Deg are degrees).
 *
 * Not covered: overlapping joints (§4.4), grouted joints (§4.5), ring-stiffened joints, the ISO 19902 §14
 * strength-factor set, thickened-can short-chord reduction (§4.3.5).
 */

function pos(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}

/**
 * Non-dimensional joint geometry and validity range of API RP 2A-WSD §4.3.1
 * (0.2 ≤ β ≤ 1.0, 10 ≤ γ ≤ 50, 30° ≤ θ ≤ 90°, Fy ≤ 500 MPa, g/D > −0.6 for K joints).
 * @param {object} o
 * @param {number} o.D chord outer diameter [m]
 * @param {number} o.T chord wall thickness [m]
 * @param {number} o.d brace outer diameter [m]
 * @param {number} o.t brace wall thickness [m]
 * @param {number} o.thetaDeg included angle between brace and chord [deg]
 * @param {number} [o.gap] gap between braces along the chord surface [m] (negative for overlap)
 * @param {number} [o.fyChord] chord yield strength [Pa]
 * @returns {{beta:number, gamma:number, tau:number, theta:number, gapRatio:(number|undefined), warnings:string[]}} β, γ, τ [-], θ [rad], g/D [-]
 */
export function jointGeometry({ D, T, d, t, thetaDeg, gap, fyChord }) {
  pos('D', D); pos('T', T); pos('d', d); pos('t', t); pos('thetaDeg', thetaDeg);
  if (d > D * (1 + 1e-9)) throw new RangeError('d must not exceed D');
  if (thetaDeg > 90) throw new RangeError('thetaDeg must not exceed 90');
  const beta = d / D; const gamma = D / (2 * T); const tau = t / T;
  const warnings = [];
  if (beta < 0.2) warnings.push('β is below 0.2 (validity range 0.2–1.0).');
  if (gamma < 10 || gamma > 50) warnings.push('γ is outside 10–50.');
  if (thetaDeg < 30) warnings.push('θ is below 30°.');
  if (fyChord !== undefined && fyChord > 500e6) warnings.push('Chord yield strength exceeds 500 MPa.');
  if (gap !== undefined && gap / D <= -0.6) warnings.push('g/D is not greater than −0.6.');
  return { beta, gamma, tau, theta: thetaDeg * Math.PI / 180, gapRatio: gap === undefined ? undefined : gap / D, warnings };
}

/**
 * Geometric factor Qβ (API RP 2A-WSD Table 4.3-1 notes): Qβ = 0.3/[β(1 − 0.833β)] for β > 0.6, else 1.0.
 * @param {number} beta diameter ratio [-]
 * @returns {number} Qβ [-]
 */
export function qBeta(beta) {
  pos('beta', beta);
  return beta > 0.6 ? 0.3 / (beta * (1 - 0.833 * beta)) : 1.0;
}

/**
 * Gap factor Qg (API RP 2A-WSD Table 4.3-1 notes): Qg = 1 + 0.2[1 − 2.8 g/D]³ ≥ 1.0 for g/D ≥ 0.05;
 * Qg = 0.13 + 0.65 φ γ^0.5 for g/D ≤ −0.05 with φ = t Fyb/(T Fy); linear interpolation between ±0.05.
 * @param {object} o
 * @param {number} o.gapRatio g/D [-]
 * @param {number} o.gamma γ [-]
 * @param {number} [o.tau] τ = t/T [-] (needed for g/D < 0.05)
 * @param {number} [o.fyBrace] brace yield strength [Pa] (needed for g/D < 0.05)
 * @param {number} [o.fyChord] chord yield strength [Pa] (needed for g/D < 0.05)
 * @returns {number} Qg [-]
 */
export function qGap({ gapRatio, gamma, tau, fyBrace, fyChord }) {
  if (!Number.isFinite(gapRatio)) throw new RangeError('gapRatio must be finite');
  pos('gamma', gamma);
  const pos05 = Math.max(1 + 0.2 * (1 - 2.8 * 0.05) ** 3, 1.0);
  if (gapRatio >= 0.05) return Math.max(1 + 0.2 * (1 - 2.8 * gapRatio) ** 3, 1.0);
  pos('tau', tau); pos('fyBrace', fyBrace); pos('fyChord', fyChord);
  const neg05 = 0.13 + 0.65 * (tau * fyBrace / fyChord) * Math.sqrt(gamma);
  if (gapRatio <= -0.05) return neg05;
  return neg05 + (pos05 - neg05) * (gapRatio + 0.05) / 0.1;
}

/**
 * Strength factor Qu (API RP 2A-WSD 21st ed. Supp. 2 / 22nd ed. Table 4.3-1).
 *  K: axial (16 + 1.2γ) β^1.2 Qg ≤ 40 β^1.2 Qg; T/Y: tension 30β, compression 2.8 + (20 + 0.8γ) β^1.6 ≤ 2.8 + 36 β^1.6;
 *  X: tension 23β (β ≤ 0.9), 20.7 + (β − 0.9)(17γ − 220) (β > 0.9), compression [2.8 + (12 + 0.1γ)β] Qβ;
 *  all types: in-plane bending (5 + 0.7γ) β^1.2, out-of-plane bending 2.5 + (4.5 + 0.2γ) β^2.6.
 * @param {object} o
 * @param {'K'|'Y'|'X'} o.type joint class for the load component ('T' is accepted as 'Y')
 * @param {'tension'|'compression'|'ipb'|'opb'} o.load brace load component
 * @param {number} o.beta β [-]
 * @param {number} o.gamma γ [-]
 * @param {number} [o.Qg=1] gap factor for K joints [-]
 * @returns {number} Qu [-]
 */
export function quApi({ type, load, beta, gamma, Qg = 1 }) {
  pos('beta', beta); pos('gamma', gamma);
  if (load === 'ipb') return (5 + 0.7 * gamma) * beta ** 1.2;
  if (load === 'opb') return 2.5 + (4.5 + 0.2 * gamma) * beta ** 2.6;
  if (load !== 'tension' && load !== 'compression') throw new RangeError(`load "${load}" is not recognised`);
  const ty = type === 'T' ? 'Y' : type;
  if (ty === 'K') return Math.min((16 + 1.2 * gamma) * beta ** 1.2 * Qg, 40 * beta ** 1.2 * Qg);
  if (ty === 'Y') {
    if (load === 'tension') return 30 * beta;
    return Math.min(2.8 + (20 + 0.8 * gamma) * beta ** 1.6, 2.8 + 36 * beta ** 1.6);
  }
  if (ty === 'X') {
    if (load === 'tension') return beta <= 0.9 ? 23 * beta : 20.7 + (beta - 0.9) * (17 * gamma - 220);
    return (2.8 + (12 + 0.1 * gamma) * beta) * qBeta(beta);
  }
  throw new RangeError(`type "${type}" is not recognised`);
}

/**
 * Chord load factor Qf = 1 + C1 (FS Pc/Py) − C2 (FS Mipb/Mp) − C3 A², A = √[(FS Pc/Py)² + (FS Mc/Mp)²]
 * (API RP 2A-WSD Eq. 4.3-2, 4.3-3; Table 4.3-2: K joints under brace axial load C = 0.2, 0.2, 0.3; T/Y joints under
 * brace axial load 0.3, 0, 0.8; X joints under brace axial load 0.2, 0, 0.5 for β ≤ 0.9 and −0.2, 0, 0.2 for β = 1.0 with
 * linear interpolation between; all joints under brace moment loading 0.2, 0, 0.4). Pc is the chord axial load
 * (tension positive), Mc = √(Mipb² + Mopb²) the chord bending resultant; the average of the chord loads on either side
 * of the brace intersection should be used. Qf is not taken above 1.0 … Qf may exceed 1.0 for chord tension per the
 * equation; the value from the equation is returned unmodified.
 * @param {object} o
 * @param {'K'|'Y'|'X'} o.type joint class
 * @param {'axial'|'moment'} o.load brace load kind
 * @param {number} o.beta β [-]
 * @param {number} o.Pc chord axial force [N], tension positive
 * @param {number} [o.Mipb=0] chord in-plane bending moment [N·m]; positive when it produces compression on the joint footprint
 * @param {number} [o.Mopb=0] chord out-of-plane bending moment [N·m]
 * @param {number} o.Py chord axial yield capacity A·Fy [N]
 * @param {number} o.Mp chord plastic moment capacity [N·m]
 * @param {number} [o.FS=1.6] safety factor [-] (1.6 normal, 1.2 where the one-third increase applies)
 * @returns {{Qf:number, A:number, C1:number, C2:number, C3:number}} factor and coefficients [-]
 */
export function qfApi({ type, load, beta, Pc, Mipb = 0, Mopb = 0, Py, Mp, FS = 1.6 }) {
  pos('Py', Py); pos('Mp', Mp); pos('FS', FS);
  let C1; let C2; let C3;
  if (load === 'moment') { C1 = 0.2; C2 = 0; C3 = 0.4; }
  else if (load === 'axial') {
    const ty = type === 'T' ? 'Y' : type;
    if (ty === 'K') { C1 = 0.2; C2 = 0.2; C3 = 0.3; }
    else if (ty === 'Y') { C1 = 0.3; C2 = 0; C3 = 0.8; }
    else if (ty === 'X') {
      const w = Math.min(Math.max((beta - 0.9) / 0.1, 0), 1);
      C1 = 0.2 + w * (-0.2 - 0.2); C2 = 0; C3 = 0.5 + w * (0.2 - 0.5);
    } else throw new RangeError(`type "${type}" is not recognised`);
  } else throw new RangeError(`load "${load}" is not recognised`);
  const pr = FS * Pc / Py;
  const mr = FS * Math.hypot(Mipb, Mopb) / Mp;
  const A = Math.hypot(pr, mr);
  return { Qf: 1 + C1 * pr - C2 * (FS * Mipb / Mp) - C3 * A * A, A, C1, C2, C3 };
}

/**
 * Classify the axial load of each brace at a planar joint into K, X and Y fractions from the load path
 * (API RP 2A-WSD §4.2.1): the component of a brace force normal to the chord is K-action as far as it is
 * balanced by braces on the same side of the chord, X-action where it is carried through the chord to braces on
 * the opposite side, and Y-action where it is reacted as beam shear in the chord. The balance is allocated in
 * proportion to the available opposing normal forces.
 * @param {Array<{P:number, thetaDeg:number, side:(1|-1)}>} braces axial force [N] (tension positive), angle to the chord [deg], side of the chord
 * @returns {Array<{K:number, X:number, Y:number}>} fractions of each brace's axial load (sum to 1)
 */
export function classifyJoint(braces) {
  if (!Array.isArray(braces) || braces.length === 0) throw new RangeError('braces must be a non-empty array');
  // normal force on the chord wall: positive = pushing into the chord (brace compression)
  const n = braces.map((b, i) => {
    if (!Number.isFinite(b.P) || !(b.thetaDeg > 0 && b.thetaDeg <= 90) || (b.side !== 1 && b.side !== -1)) {
      throw new RangeError(`braces[${i}] needs finite P, 0 < thetaDeg <= 90 and side = ±1`);
    }
    return -b.P * Math.sin(b.thetaDeg * Math.PI / 180);
  });
  const out = braces.map(() => ({ K: 0, X: 0, Y: 1 }));
  for (const side of [1, -1]) {
    let push = 0; let pull = 0;
    braces.forEach((b, i) => { if (b.side === side) { if (n[i] > 0) push += n[i]; else pull -= n[i]; } });
    const bal = Math.min(push, pull);
    braces.forEach((b, i) => {
      if (b.side !== side || n[i] === 0) return;
      const tot = n[i] > 0 ? push : pull;
      out[i].K = tot > 0 ? bal / tot : 0;
    });
  }
  // residual normal force per side (after K balance), acting on the chord in the +side direction sense
  const resid = { 1: 0, '-1': 0 };
  braces.forEach((b, i) => { resid[b.side] += n[i]; });
  // X-action: residual push on one side balanced by residual push on the other (crushing) or pull by pull
  const r1 = resid[1]; const r2 = resid[-1];
  let xbal = 0;
  if (r1 * r2 > 0) xbal = Math.min(Math.abs(r1), Math.abs(r2));
  braces.forEach((b, i) => {
    const r = resid[b.side];
    if (n[i] === 0) { out[i] = { K: 0, X: 0, Y: 1 }; return; }
    const rest = 1 - out[i].K;
    const sameSense = n[i] * r > 0;
    const xfrac = sameSense && Math.abs(r) > 0 ? rest * xbal / Math.abs(r) : 0;
    out[i].X = xfrac;
    out[i].Y = Math.max(0, 1 - out[i].K - xfrac);
  });
  return out;
}

// Utilisation reported when the chord load factor leaves the joint without capacity.
const OVERSTRESSED_UC = 1000;

/**
 * Allowable capacities and interaction check of a simple tubular joint (API RP 2A-WSD §4.3):
 * Pa = Qu Qf Fyc T²/(FS sinθ), Ma = Qu Qf Fyc T² d/(FS sinθ) (Eq. 4.3-1a/b); interaction
 * |P/Pa| + (M/Ma)²ipb + |M/Ma|opb ≤ 1.0 (Eq. 4.3-5); minimum capacity: the joint should develop at least 50 % of the
 * effective strength of the brace (§4.2.3). Mixed classification uses the weighted average of the K, X and Y axial
 * capacities (§4.3.2). Fyc is limited to 0.8 of the chord tensile strength when `fuChord` is given.
 * Overstressed chord: Qf = 1 + C1 (FS Pc/Py) − C2 (FS Mipb/Mp) − C3 A² turns zero and then negative when the
 * chord's own utilisation A is high enough. A negative factor would give negative capacities and, through
 * |P|/Pa, a negative — apparently safe — utilisation. The joint has no capacity left in that state: Qf is
 * floored at zero, the affected capacities are zero, the corresponding utilisation terms and the interaction
 * value are set to 1000 (a large finite number that ranks and plots), and `chordOverstressed` is true.
 * @param {object} o
 * @param {number} o.D chord outer diameter [m]
 * @param {number} o.T chord wall thickness [m]
 * @param {number} o.d brace outer diameter [m]
 * @param {number} o.t brace wall thickness [m]
 * @param {number} o.thetaDeg brace angle [deg]
 * @param {number} o.fyChord chord yield strength [Pa]
 * @param {number} [o.fuChord] chord tensile strength [Pa]
 * @param {number} [o.fyBrace=o.fyChord] brace yield strength [Pa]
 * @param {{K?:number, X?:number, Y?:number}} [o.classification={Y:1}] fractions of the brace axial load in each class
 * @param {number} [o.gap] gap for the K portion [m]
 * @param {number} [o.P=0] brace axial force [N], tension positive
 * @param {number} [o.Mipb=0] brace in-plane bending moment [N·m]
 * @param {number} [o.Mopb=0] brace out-of-plane bending moment [N·m]
 * @param {{P?:number, Mipb?:number, Mopb?:number}} [o.chord] chord loads at the joint [N, N·m] (P tension positive)
 * @param {number} [o.FS=1.6] safety factor [-]
 * @param {number} [o.braceStrength] effective brace strength for the 50 % rule [N]; default A_brace·Fyb
 * @returns {{geometry:object, Qu:object, Qf:object, Pa:number, MaIpb:number, MaOpb:number, uc:number, ucAxial:number, ucIpb:number, ucOpb:number,
 *   minimumCapacity:{required:number, provided:number, ok:boolean}, checks:Array, warnings:string[], chordOverstressed:boolean}}
 *   capacities [N, N·m] (never negative), utilisations [-] (finite); Qf as floored at zero; chordOverstressed = a
 *   chord load factor was not positive
 */
export function apiJointCheck({
  D, T, d, t, thetaDeg, fyChord, fuChord, fyBrace, classification = { Y: 1 }, gap, P = 0, Mipb = 0, Mopb = 0, chord = {},
  FS = 1.6, braceStrength,
}) {
  pos('fyChord', fyChord);
  const fyb = fyBrace === undefined ? fyChord : pos('fyBrace', fyBrace);
  const geo = jointGeometry({ D, T, d, t, thetaDeg, gap, fyChord });
  const warnings = [...geo.warnings];
  const Fyc = fuChord !== undefined ? Math.min(fyChord, 0.8 * fuChord) : fyChord;
  const { beta, gamma, tau, theta } = geo;
  const fK = classification.K || 0; const fX = classification.X || 0; const fY = classification.Y || 0;
  const fsum = fK + fX + fY;
  if (!(fsum > 0)) throw new RangeError('classification fractions must sum to a positive number');
  const di = D - 2 * T;
  const Ac = Math.PI / 4 * (D * D - di * di);
  const Py = Ac * Fyc;
  const Mp = (D ** 3 - di ** 3) / 6 * Fyc;
  const load = P >= 0 ? 'tension' : 'compression';
  const base = Fyc * T * T / (FS * Math.sin(theta));
  const Qu = {}; const Qf = {};
  let Pa = 0;
  let chordOverstressed = false;
  for (const [ty, fr] of [['K', fK], ['X', fX], ['Y', fY]]) {
    if (!(fr > 0)) continue;
    let Qg = 1;
    if (ty === 'K') {
      if (gap === undefined) throw new RangeError('gap is required for the K portion of the joint');
      Qg = qGap({ gapRatio: gap / D, gamma, tau, fyBrace: fyb, fyChord: Fyc });
    }
    Qu[ty] = quApi({ type: ty, load, beta, gamma, Qg });
    Qf[ty] = Math.max(qfApi({ type: ty, load: 'axial', beta, Pc: chord.P || 0, Mipb: chord.Mipb || 0, Mopb: chord.Mopb || 0, Py, Mp, FS }).Qf, 0);
    if (!(Qf[ty] > 0)) chordOverstressed = true;
    Pa += (fr / fsum) * Qu[ty] * Qf[ty] * base;
  }
  Qu.ipb = quApi({ type: 'Y', load: 'ipb', beta, gamma });
  Qu.opb = quApi({ type: 'Y', load: 'opb', beta, gamma });
  Qf.moment = Math.max(qfApi({ type: 'Y', load: 'moment', beta, Pc: chord.P || 0, Mipb: chord.Mipb || 0, Mopb: chord.Mopb || 0, Py, Mp, FS }).Qf, 0);
  if (!(Qf.moment > 0)) chordOverstressed = true;
  const MaIpb = Qu.ipb * Qf.moment * base * d;
  const MaOpb = Qu.opb * Qf.moment * base * d;
  if (chordOverstressed) {
    warnings.push('Chord utilisation is so high that the chord load factor Qf is not positive: the joint has no capacity '
      + `and the utilisation is reported as ${OVERSTRESSED_UC}.`);
  }
  // a capacity of zero cannot carry anything: the term is the overstress marker, not a division by zero
  const term = (load, cap) => (cap > 0 ? Math.abs(load) / cap : OVERSTRESSED_UC);
  const ucAxial = term(P, Pa);
  const ucIpb = term(Mipb, MaIpb);
  const ucOpb = term(Mopb, MaOpb);
  const uc = chordOverstressed ? OVERSTRESSED_UC : ucAxial + ucIpb * ucIpb + ucOpb;
  const dbi = d - 2 * t;
  const braceCap = braceStrength === undefined ? Math.PI / 4 * (d * d - dbi * dbi) * fyb : pos('braceStrength', braceStrength);
  // the 50 % rule compares ultimate values: allowable × FS against half of the brace strength
  const provided = Pa * FS;
  const checks = [
    { name: 'Joint interaction', uc, ref: 'API RP 2A-WSD Eq. 4.3-5' },
    { name: 'Minimum joint capacity (50 % of brace strength)', uc: provided > 0 ? 0.5 * braceCap / provided : OVERSTRESSED_UC, ref: 'API RP 2A-WSD §4.2.3' },
  ];
  return {
    geometry: geo, Qu, Qf, Pa, MaIpb, MaOpb, uc, ucAxial, ucIpb, ucOpb,
    minimumCapacity: { required: 0.5 * braceCap, provided, ok: provided >= 0.5 * braceCap }, checks, warnings, chordOverstressed,
  };
}
