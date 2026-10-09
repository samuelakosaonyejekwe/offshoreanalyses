/**
 * Stress concentration factors for simple tubular joints (Efthymiou equations as tabulated in
 * DNV-RP-C203 Appendix B, Tables B-1 to B-5) and for girth welds (DNV-RP-C203 §3.3.7), plus the superposition
 * of hot-spot stresses at eight points around the brace–chord intersection (DNV-RP-C203 §3.3.2).
 *
 * Parameters: β = d/D, γ = D/(2T), τ = t/T, α = 2L/D, ζ = g/D, θ = brace angle. Equation numbers quoted below
 * are those of DNV-RP-C203 Appendix B. Validity range (enforced with warnings, values are still returned):
 * 0.2 ≤ β ≤ 1.0, 0.2 ≤ τ ≤ 1.0, 8 ≤ γ ≤ 32, 4 ≤ α ≤ 40, 20° ≤ θ ≤ 90°, −0.6β/sinθ ≤ ζ ≤ 1.0.
 *
 * Not covered: KT joints with a single brace loaded, tubular-to-cone transitions, ring-stiffened joints,
 * the minimum SCF and effective-SCF rules for saddle-dominated in-plane bending.
 */

function pos(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}

/**
 * Non-dimensional parameters of a brace–chord intersection with the Efthymiou validity check.
 * @param {object} o
 * @param {number} o.D chord outer diameter [m]
 * @param {number} o.T chord wall thickness [m]
 * @param {number} o.d brace outer diameter [m]
 * @param {number} o.t brace wall thickness [m]
 * @param {number} o.L chord length between supports / points of contraflexure [m]
 * @param {number} o.thetaDeg brace angle [deg]
 * @param {number} [o.gap] gap to the neighbouring brace [m]
 * @returns {{beta:number, gamma:number, tau:number, alpha:number, theta:number, zeta:(number|undefined), warnings:string[]}}
 */
export function scfParameters({ D, T, d, t, L, thetaDeg, gap }) {
  pos('D', D); pos('T', T); pos('d', d); pos('t', t); pos('L', L); pos('thetaDeg', thetaDeg);
  if (thetaDeg > 90) throw new RangeError('thetaDeg must not exceed 90');
  const p = { beta: d / D, gamma: D / (2 * T), tau: t / T, alpha: 2 * L / D, theta: thetaDeg * Math.PI / 180 };
  p.zeta = gap === undefined ? undefined : gap / D;
  p.warnings = validity(p);
  return p;
}

function validity({ beta, gamma, tau, alpha, theta, zeta }) {
  const w = [];
  if (beta < 0.2 || beta > 1.0) w.push(`β = ${beta.toFixed(3)} is outside 0.2–1.0.`);
  if (tau < 0.2 || tau > 1.0) w.push(`τ = ${tau.toFixed(3)} is outside 0.2–1.0.`);
  if (gamma < 8 || gamma > 32) w.push(`γ = ${gamma.toFixed(2)} is outside 8–32.`);
  if (alpha !== undefined && (alpha < 4 || alpha > 40)) w.push(`α = ${alpha.toFixed(2)} is outside 4–40.`);
  const deg = theta * 180 / Math.PI;
  if (deg < 20 - 1e-9 || deg > 90 + 1e-9) w.push(`θ = ${deg.toFixed(1)}° is outside 20°–90°.`);
  if (zeta !== undefined && (zeta > 1.0 || zeta < -0.6 * beta / Math.sin(theta))) w.push(`ζ = ${zeta.toFixed(3)} is outside −0.6β/sinθ … 1.0.`);
  return w;
}

/**
 * Short-chord correction factors (α < 12), DNV-RP-C203 Table B-1:
 * F1 = 1 − (0.83β − 0.56β² − 0.02) γ^0.23 exp(−0.21 γ^−1.16 α^2.5),
 * F2 = 1 − (1.43β − 0.97β² − 0.03) γ^0.04 exp(−0.71 γ^−1.38 α^2.5),
 * F3 = 1 − 0.55 β^1.8 γ^0.16 exp(−0.49 γ^−0.89 α^1.8),
 * F4 = 1 − 1.07 β^1.88 exp(−0.16 γ^−1.06 α^2.4) (Table B-3). All equal 1.0 for α ≥ 12.
 * @param {object} o
 * @param {number} o.beta β [-]
 * @param {number} o.gamma γ [-]
 * @param {number} o.alpha α [-]
 * @returns {{F1:number, F2:number, F3:number, F4:number}}
 */
export function shortChordFactors({ beta, gamma, alpha }) {
  pos('beta', beta); pos('gamma', gamma); pos('alpha', alpha);
  if (alpha >= 12) return { F1: 1, F2: 1, F3: 1, F4: 1 };
  return {
    F1: 1 - (0.83 * beta - 0.56 * beta * beta - 0.02) * gamma ** 0.23 * Math.exp(-0.21 * gamma ** -1.16 * alpha ** 2.5),
    F2: 1 - (1.43 * beta - 0.97 * beta * beta - 0.03) * gamma ** 0.04 * Math.exp(-0.71 * gamma ** -1.38 * alpha ** 2.5),
    F3: 1 - 0.55 * beta ** 1.8 * gamma ** 0.16 * Math.exp(-0.49 * gamma ** -0.89 * alpha ** 1.8),
    F4: 1 - 1.07 * beta ** 1.88 * Math.exp(-0.16 * gamma ** -1.06 * alpha ** 2.4),
  };
}

/* Shared building blocks (equation numbers of DNV-RP-C203 Appendix B). */
const eq8 = (b, g, t, th) => 1.45 * b * t ** 0.85 * g ** (1 - 0.68 * b) * Math.sin(th) ** 0.7;
const eq9 = (b, g, t, th) => 1 + 0.65 * b * t ** 0.4 * g ** (1.09 - 0.77 * b) * Math.sin(th) ** (0.06 * g - 1.16);
const eq10 = (b, g, t, th) => g * t * b * (1.7 - 1.05 * b ** 3) * Math.sin(th) ** 1.6;
const braceFromChordOpb = (b, g, t) => t ** -0.54 * g ** -0.05 * (0.99 - 0.47 * b + 0.08 * b ** 4);

/**
 * Efthymiou SCFs for a simple T or Y joint (DNV-RP-C203 Table B-1).
 * Axial load, chord ends fixed: Eq. (1)–(4); general fixity: Eq. (5), (6), (3), (7) with C1 = 2(C − 0.5), C2 = C/2,
 * C3 = C/5 and 0.5 ≤ C ≤ 1.0 (typically 0.7). In-plane bending: Eq. (8), (9). Out-of-plane bending: Eq. (10), (11).
 * Short-chord factors: F1 on the saddle SCFs for fixed ends, F2 for general fixity, F3 for out-of-plane bending.
 * @param {object} o
 * @param {number} o.beta β [-]
 * @param {number} o.gamma γ [-]
 * @param {number} o.tau τ [-]
 * @param {number} o.alpha α = 2L/D [-]
 * @param {number} o.theta brace angle [rad]
 * @param {'fixed'|'general'} [o.chordEnds='general'] chord end fixity model for axial load
 * @param {number} [o.C=0.7] chord end fixity parameter [-], 0.5 (pinned) … 1.0 (fixed), general fixity only
 * @returns {{axial:{chordSaddle:number, chordCrown:number, braceSaddle:number, braceCrown:number}, ipb:{chordCrown:number, braceCrown:number},
 *   opb:{chordSaddle:number, braceSaddle:number}, factors:object, warnings:string[]}} SCFs [-]
 */
export function scfTY({ beta, gamma, tau, alpha, theta, chordEnds = 'general', C = 0.7 }) {
  pos('beta', beta); pos('gamma', gamma); pos('tau', tau); pos('alpha', alpha); pos('theta', theta);
  if (chordEnds !== 'fixed' && chordEnds !== 'general') throw new RangeError('chordEnds must be "fixed" or "general"');
  if (!(C >= 0.5 && C <= 1.0)) throw new RangeError('C must lie in [0.5, 1.0]');
  const s = Math.sin(theta);
  const F = shortChordFactors({ beta, gamma, alpha });
  const e1 = gamma * tau ** 1.1 * (1.11 - 3 * (beta - 0.52) ** 2) * s ** 1.6;
  const e3 = 1.3 + gamma * tau ** 0.52 * alpha ** 0.1 * (0.187 - 1.25 * beta ** 1.1 * (beta - 0.96)) * s ** (2.7 - 0.01 * alpha);
  const crownBase = gamma ** 0.2 * tau * (2.65 + 5 * (beta - 0.65) ** 2);
  const braceCrownBase = 3 + gamma ** 1.2 * (0.12 * Math.exp(-4 * beta) + 0.011 * beta * beta - 0.045);
  let axial;
  if (chordEnds === 'fixed') {
    axial = {
      chordSaddle: e1 * F.F1,
      chordCrown: crownBase + tau * beta * (0.25 * alpha - 3) * s,
      braceSaddle: e3 * F.F1,
      braceCrown: braceCrownBase + beta * tau * (0.1 * alpha - 1.2),
    };
  } else {
    const C1 = 2 * (C - 0.5); const C2 = C / 2; const C3 = C / 5;
    axial = {
      chordSaddle: (e1 + C1 * (0.8 * alpha - 6) * tau * beta * beta * Math.sqrt(1 - beta * beta) * Math.sin(2 * theta) ** 2) * F.F2,
      chordCrown: crownBase + tau * beta * (C2 * alpha - 3) * s,
      braceSaddle: e3 * F.F2,
      braceCrown: braceCrownBase + beta * tau * (C3 * alpha - 1.2),
    };
  }
  const e10 = eq10(beta, gamma, tau, theta);
  return {
    axial,
    ipb: { chordCrown: eq8(beta, gamma, tau, theta), braceCrown: eq9(beta, gamma, tau, theta) },
    opb: { chordSaddle: e10 * F.F3, braceSaddle: braceFromChordOpb(beta, gamma, tau) * e10 * F.F3 },
    factors: F, warnings: validity({ beta, gamma, tau, alpha, theta }),
  };
}

/**
 * Efthymiou SCFs for an X joint (DNV-RP-C203 Table B-2). Balanced axial load: Eq. (12)–(15); in-plane bending
 * Eq. (8), (9); balanced out-of-plane bending Eq. (16), (17). Short chords (α < 12): F1 (chord ends fixed) or F2
 * (pinned) on the axial saddle SCFs, F3 on the out-of-plane SCFs.
 * @param {object} o β, γ, τ, α, θ [rad] as scfTY
 * @param {'fixed'|'pinned'} [o.chordEnds='fixed'] selects F1 or F2 for short chords
 * @returns {{axial:object, ipb:object, opb:object, factors:object, warnings:string[]}} SCFs [-]
 */
export function scfX({ beta, gamma, tau, alpha, theta, chordEnds = 'fixed' }) {
  pos('beta', beta); pos('gamma', gamma); pos('tau', tau); pos('alpha', alpha); pos('theta', theta);
  const s = Math.sin(theta);
  const F = shortChordFactors({ beta, gamma, alpha });
  const Fs = chordEnds === 'fixed' ? F.F1 : F.F2;
  const e16 = gamma * tau * beta * (1.56 - 1.34 * beta ** 4) * s ** 1.6;
  return {
    axial: {
      chordSaddle: 3.87 * gamma * tau * beta * (1.10 - beta ** 1.8) * s ** 1.7 * Fs,
      chordCrown: gamma ** 0.2 * tau * (2.65 + 5 * (beta - 0.65) ** 2) - 3 * tau * beta * s,
      braceSaddle: (1 + 1.9 * gamma * tau ** 0.5 * beta ** 0.9 * (1.09 - beta ** 1.7) * s ** 2.5) * Fs,
      braceCrown: 3 + gamma ** 1.2 * (0.12 * Math.exp(-4 * beta) + 0.011 * beta * beta - 0.045),
    },
    ipb: { chordCrown: eq8(beta, gamma, tau, theta), braceCrown: eq9(beta, gamma, tau, theta) },
    opb: { chordSaddle: e16 * F.F3, braceSaddle: braceFromChordOpb(beta, gamma, tau) * e16 * F.F3 },
    factors: F, warnings: validity({ beta, gamma, tau, alpha, theta }),
  };
}

/**
 * Efthymiou SCFs for a gap or overlap K joint (DNV-RP-C203 Table B-3) for the brace under consideration (A)
 * with neighbouring brace B.
 * Balanced axial load: chord Eq. (20)
 *   τ^0.9 γ^0.5 (0.67 − β² + 1.16β) sinθ (sinθmax/sinθmin)^0.30 (βmax/βmin)^0.30 [1.64 + 0.29 β^−0.38 atan(8ζ)],
 * brace Eq. (21)
 *   1 + (1.97 − 1.57 β^0.25) τ^−0.14 (sinθ)^0.7 · Eq.(20) + sin^1.8(θA + θB) [0.131 − 0.084 atan(14ζ + 4.2β)] C β^1.5 γ^0.5 τ^−1.22
 * with C = 0 for gap joints, 1 for the through brace and 0.5 for the overlapping brace.
 * Unbalanced in-plane bending: Eq. (8), (9). Unbalanced out-of-plane bending: chord saddle Eq. (22)
 *   [Eq.10]A [1 − 0.08 (βB γ)^0.5 exp(−0.8x)] + [Eq.10]B [1 − 0.08 (βA γ)^0.5 exp(−0.8x)] [2.05 βmax^0.5 exp(−1.3x)], x = 1 + ζ sinθA/βA,
 * brace saddle Eq. (23) = τ^−0.54 γ^−0.05 (0.99 − 0.47β + 0.08β⁴) · Eq.(22); F4 applies to Eq. (10) for short chords.
 * @param {object} o
 * @param {number} o.gamma γ [-]
 * @param {number} o.alpha α [-]
 * @param {number} o.zeta ζ = g/D [-] (negative for overlap)
 * @param {{beta:number, tau:number, theta:number}} o.A brace under consideration (θ in rad)
 * @param {{beta:number, tau:number, theta:number}} o.B the other brace
 * @param {number} [o.C=0] 0 gap joint, 1 through brace, 0.5 overlapping brace
 * @returns {{axial:{chord:number, brace:number}, ipb:{chordCrown:number, braceCrown:number}, opb:{chordSaddle:number, braceSaddle:number},
 *   warnings:string[]}} SCFs [-] for brace A (axial SCFs apply at the governing location, crown or saddle)
 */
export function scfK({ gamma, alpha, zeta, A, B, C = 0 }) {
  pos('gamma', gamma); pos('alpha', alpha);
  if (!Number.isFinite(zeta)) throw new RangeError('zeta must be finite');
  for (const [nm, br] of [['A', A], ['B', B]]) {
    if (!br) throw new RangeError(`${nm} must be given`);
    pos(`${nm}.beta`, br.beta); pos(`${nm}.tau`, br.tau); pos(`${nm}.theta`, br.theta);
  }
  if (![0, 0.5, 1].includes(C)) throw new RangeError('C must be 0, 0.5 or 1');
  const sA = Math.sin(A.theta); const sB = Math.sin(B.theta);
  const smax = Math.max(sA, sB); const smin = Math.min(sA, sB);
  const bmax = Math.max(A.beta, B.beta); const bmin = Math.min(A.beta, B.beta);
  const b = A.beta; const t = A.tau;
  const e20 = t ** 0.9 * Math.sqrt(gamma) * (0.67 - b * b + 1.16 * b) * sA * (smax / smin) ** 0.30 * (bmax / bmin) ** 0.30
    * (1.64 + 0.29 * b ** -0.38 * Math.atan(8 * zeta));
  const e21 = 1 + (1.97 - 1.57 * b ** 0.25) * t ** -0.14 * sA ** 0.7 * e20
    + Math.sin(A.theta + B.theta) ** 1.8 * (0.131 - 0.084 * Math.atan(14 * zeta + 4.2 * b)) * C * b ** 1.5 * Math.sqrt(gamma) * t ** -1.22;
  const F4a = shortChordFactors({ beta: A.beta, gamma, alpha }).F4;
  const F4b = shortChordFactors({ beta: B.beta, gamma, alpha }).F4;
  const x = 1 + zeta * sA / A.beta;
  const e10A = eq10(A.beta, gamma, A.tau, A.theta) * F4a;
  const e10B = eq10(B.beta, gamma, B.tau, B.theta) * F4b;
  const e22 = e10A * (1 - 0.08 * Math.sqrt(B.beta * gamma) * Math.exp(-0.8 * x))
    + e10B * (1 - 0.08 * Math.sqrt(A.beta * gamma) * Math.exp(-0.8 * x)) * (2.05 * Math.sqrt(bmax) * Math.exp(-1.3 * x));
  const warnings = [
    ...validity({ beta: A.beta, gamma, tau: A.tau, alpha, theta: A.theta, zeta }),
    ...validity({ beta: B.beta, gamma, tau: B.tau, alpha, theta: B.theta, zeta }),
  ];
  return {
    axial: { chord: e20, brace: e21 },
    ipb: { chordCrown: eq8(b, gamma, t, A.theta), braceCrown: eq9(b, gamma, t, A.theta) },
    opb: { chordSaddle: e22, braceSaddle: braceFromChordOpb(b, gamma, t) * e22 },
    warnings: [...new Set(warnings)],
  };
}

/**
 * Efthymiou SCFs for a KT joint with all braces loaded (DNV-RP-C203 Table B-4): diagonal braces A and C, central
 * brace B.
 * Balanced axial load: Eq. (20) and (21) of the K joint with, for the diagonal braces, ζ = ζAB + ζBC + βB and, for
 * the central brace, ζ = max(ζAB, ζBC); the neighbouring brace for the βmax/βmin and θ terms is the other diagonal
 * (for A, C) or the diagonal across the larger gap (for B).
 * Unbalanced out-of-plane bending, chord saddle adjacent to diagonal brace A:
 *   [Eq.10]A [1 − 0.08(βBγ)^0.5 e^(−0.8xAB)] [1 − 0.08(βCγ)^0.5 e^(−0.8xAC)] + [Eq.10]B [1 − 0.08(βAγ)^0.5 e^(−0.8xAB)] [2.05 βmax^0.5 e^(−1.3xAB)]
 *   + [Eq.10]C [1 − 0.08(βAγ)^0.5 e^(−0.8xAC)] [2.05 βmax^0.5 e^(−1.3xAC)], xAB = 1 + ζAB sinθA/βA, xAC = 1 + (ζAB + ζBC + βB) sinθA/βA;
 * adjacent to central brace B:
 *   [Eq.10]B [1 − 0.08(βAγ)^0.5 e^(−0.8xAB)]^((βA/βB)²) [1 − 0.08(βCγ)^0.5 e^(−0.8xBC)]^((βC/βB)²)
 *   + [Eq.10]A [1 − 0.08(βBγ)^0.5 e^(−0.8xAB)] [2.05 βmax^0.5 e^(−1.3xAB)] + [Eq.10]C [1 − 0.08(βBγ)^0.5 e^(−0.8xBC)] [2.05 βmax^0.5 e^(−1.3xBC)],
 *   xAB = 1 + ζAB sinθB/βB, xBC = 1 + ζBC sinθB/βB;
 * brace saddle = τ^−0.54 γ^−0.05 (0.99 − 0.47β + 0.08β⁴) × chord saddle. In-plane bending: Eq. (8), (9).
 * @param {object} o
 * @param {number} o.gamma γ [-]
 * @param {number} o.alpha α [-]
 * @param {number} o.zetaAB gap ratio between A and B [-]
 * @param {number} o.zetaBC gap ratio between B and C [-]
 * @param {{beta:number, tau:number, theta:number}} o.A diagonal brace
 * @param {{beta:number, tau:number, theta:number}} o.B central brace
 * @param {{beta:number, tau:number, theta:number}} o.C diagonal brace
 * @returns {{A:object, B:object, C:object, warnings:string[]}} per brace: { axial:{chord, brace}, ipb:{chordCrown, braceCrown}, opb:{chordSaddle, braceSaddle} }
 */
export function scfKT({ gamma, alpha, zetaAB, zetaBC, A, B, C }) {
  pos('gamma', gamma); pos('alpha', alpha);
  if (!Number.isFinite(zetaAB) || !Number.isFinite(zetaBC)) throw new RangeError('zetaAB and zetaBC must be finite');
  for (const [nm, br] of [['A', A], ['B', B], ['C', C]]) {
    if (!br) throw new RangeError(`${nm} must be given`);
    pos(`${nm}.beta`, br.beta); pos(`${nm}.tau`, br.tau); pos(`${nm}.theta`, br.theta);
  }
  const bmax = Math.max(A.beta, B.beta, C.beta);
  const e10 = (br) => eq10(br.beta, gamma, br.tau, br.theta) * shortChordFactors({ beta: br.beta, gamma, alpha }).F4;
  const red = (bOther, x) => 1 - 0.08 * Math.sqrt(bOther * gamma) * Math.exp(-0.8 * x);
  const carry = (x) => 2.05 * Math.sqrt(bmax) * Math.exp(-1.3 * x);
  const diag = (P, Q, zPB, zBQ) => {
    // P = diagonal under consideration, Q = the other diagonal
    const zeta = zPB + zBQ + B.beta;
    const k = scfK({ gamma, alpha, zeta, A: P, B: Q });
    const xPB = 1 + zPB * Math.sin(P.theta) / P.beta;
    const xPQ = 1 + zeta * Math.sin(P.theta) / P.beta;
    const chord = e10(P) * red(B.beta, xPB) * red(Q.beta, xPQ) + e10(B) * red(P.beta, xPB) * carry(xPB) + e10(Q) * red(P.beta, xPQ) * carry(xPQ);
    return { axial: k.axial, ipb: k.ipb, opb: { chordSaddle: chord, braceSaddle: braceFromChordOpb(P.beta, gamma, P.tau) * chord } };
  };
  const zc = Math.max(zetaAB, zetaBC);
  const kB = scfK({ gamma, alpha, zeta: zc, A: B, B: zetaAB >= zetaBC ? A : C });
  const xAB = 1 + zetaAB * Math.sin(B.theta) / B.beta;
  const xBC = 1 + zetaBC * Math.sin(B.theta) / B.beta;
  const chordB = e10(B) * red(A.beta, xAB) ** ((A.beta / B.beta) ** 2) * red(C.beta, xBC) ** ((C.beta / B.beta) ** 2)
    + e10(A) * red(B.beta, xAB) * carry(xAB) + e10(C) * red(B.beta, xBC) * carry(xBC);
  const w = new Set();
  for (const br of [A, B, C]) for (const m of validity({ beta: br.beta, gamma, tau: br.tau, alpha, theta: br.theta })) w.add(m);
  return {
    A: diag(A, C, zetaAB, zetaBC),
    B: { axial: kB.axial, ipb: kB.ipb, opb: { chordSaddle: chordB, braceSaddle: braceFromChordOpb(B.beta, gamma, B.tau) * chordB } },
    C: diag(C, A, zetaBC, zetaAB),
    warnings: [...w],
  };
}

/**
 * Hot-spot stresses at eight points around the brace–chord intersection by superposition of the axial,
 * in-plane and out-of-plane components (DNV-RP-C203 §3.3.2, Eq. 3.3.1): points 1 and 5 are the crowns, 3 and 7
 * the saddles.
 *   σ1 = SCF_AC σx + SCF_MIP σmy
 *   σ2 = ½(SCF_AC + SCF_AS) σx + ½√2 SCF_MIP σmy − ½√2 SCF_MOP σmz
 *   σ3 = SCF_AS σx − SCF_MOP σmz
 *   σ4 = ½(SCF_AC + SCF_AS) σx − ½√2 SCF_MIP σmy − ½√2 SCF_MOP σmz
 *   σ5 = SCF_AC σx − SCF_MIP σmy
 *   σ6 = ½(SCF_AC + SCF_AS) σx − ½√2 SCF_MIP σmy + ½√2 SCF_MOP σmz
 *   σ7 = SCF_AS σx + SCF_MOP σmz
 *   σ8 = ½(SCF_AC + SCF_AS) σx + ½√2 SCF_MIP σmy + ½√2 SCF_MOP σmz
 * @param {object} o
 * @param {number} o.sigmaX nominal axial stress (or stress range) in the brace [Pa]
 * @param {number} o.sigmaMy nominal in-plane bending stress [Pa]
 * @param {number} o.sigmaMz nominal out-of-plane bending stress [Pa]
 * @param {number} o.scfAC axial SCF at the crown [-]
 * @param {number} o.scfAS axial SCF at the saddle [-]
 * @param {number} o.scfMIP in-plane bending SCF [-]
 * @param {number} o.scfMOP out-of-plane bending SCF [-]
 * @returns {{sigma:number[], max:number, maxAbs:number, governingPoint:number}} hot-spot stresses at points 1–8 [Pa]; governingPoint is 1-based
 */
export function hotSpotStresses({ sigmaX = 0, sigmaMy = 0, sigmaMz = 0, scfAC, scfAS, scfMIP, scfMOP }) {
  for (const [k, v] of Object.entries({ scfAC, scfAS, scfMIP, scfMOP })) {
    if (!Number.isFinite(v)) throw new RangeError(`${k} must be a finite number`);
  }
  const h = Math.SQRT2 / 2;
  const mean = 0.5 * (scfAC + scfAS) * sigmaX;
  const ip = scfMIP * sigmaMy; const op = scfMOP * sigmaMz;
  const sigma = [
    scfAC * sigmaX + ip,
    mean + h * ip - h * op,
    scfAS * sigmaX - op,
    mean - h * ip - h * op,
    scfAC * sigmaX - ip,
    mean - h * ip + h * op,
    scfAS * sigmaX + op,
    mean + h * ip + h * op,
  ];
  let gp = 0;
  for (let i = 1; i < 8; i++) if (Math.abs(sigma[i]) > Math.abs(sigma[gp])) gp = i;
  return { sigma, max: Math.max(...sigma), maxAbs: Math.abs(sigma[gp]), governingPoint: gp + 1 };
}

/**
 * SCF for a butt weld in plates with axial misalignment (DNV-RP-C203 §3.1.2): SCF = 1 + 3 (δm − δ0)/t, δ0 = 0.1 t being
 * the misalignment already contained in the S-N curve. Not less than 1.0.
 * @param {object} o
 * @param {number} o.t plate thickness [m]
 * @param {number} o.deltaM axial misalignment (eccentricity) [m]
 * @param {number} [o.delta0=0.1·t] misalignment inherent in the S-N data [m]
 * @returns {number} SCF [-]
 */
export function scfButtWeldPlate({ t, deltaM, delta0 }) {
  pos('t', t);
  if (!(deltaM >= 0)) throw new RangeError('deltaM must be non-negative');
  const d0 = delta0 === undefined ? 0.1 * t : delta0;
  return Math.max(1, 1 + 3 * (deltaM - d0) / t);
}

/**
 * SCF at girth (butt) welds in tubulars (DNV-RP-C203 §3.3.7).
 * Equal thickness: SCF = 1 + 3 (δm − δ0)/t · exp(−√(t/D)) (set δ0 = 0 for the plain form 1 + 3 δm/t · e^(−√(t/D))).
 * Thickness transition T > t: SCF = 1 + 6 (δt + δm − δ0)/t · 1/[1 + (T/t)^β] · e^(−α), α = 1.82 L/√(D t) · 1/[1 + (T/t)^β],
 * β = 1.5 − 1.0/log10(D/t) + 3.0/[log10(D/t)]², δt = ½(T − t), L = width of the weld at the surface.
 * Not less than 1.0.
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.t (thinner) wall thickness [m]
 * @param {number} [o.T=o.t] thicker wall thickness [m]
 * @param {number} o.deltaM fabrication misalignment (out-of-roundness / centre eccentricity) [m]
 * @param {number} [o.delta0=0] misalignment inherent in the S-N curve [m] (0.1 t for curves that include it)
 * @param {number} [o.L] weld cap width [m], required for a thickness transition
 * @returns {{scf:number, alpha:(number|undefined), beta:(number|undefined), deltaT:number}} SCF [-] and intermediate factors
 */
export function scfGirthWeld({ D, t, T, deltaM, delta0 = 0, L }) {
  pos('D', D); pos('t', t);
  if (!(deltaM >= 0)) throw new RangeError('deltaM must be non-negative');
  const Tt = T === undefined ? t : pos('T', T);
  if (Tt < t) throw new RangeError('T must not be smaller than t');
  if (Tt === t) {
    return { scf: Math.max(1, 1 + 3 * (deltaM - delta0) / t * Math.exp(-Math.sqrt(t / D))), deltaT: 0 };
  }
  pos('L', L);
  const lg = Math.log10(D / t);
  const beta = 1.5 - 1.0 / lg + 3.0 / (lg * lg);
  const k = 1 / (1 + (Tt / t) ** beta);
  const alpha = 1.82 * L / Math.sqrt(D * t) * k;
  const deltaT = 0.5 * (Tt - t);
  return { scf: Math.max(1, 1 + 6 * (deltaT + deltaM - delta0) / t * k * Math.exp(-alpha)), alpha, beta, deltaT };
}

/**
 * Convenience wrapper: SCFs and the eight hot-spot stresses (chord side and brace side) of a T/Y joint from
 * dimensions and nominal brace stresses.
 * @param {object} o
 * @param {number} o.D chord diameter [m]
 * @param {number} o.T chord thickness [m]
 * @param {number} o.d brace diameter [m]
 * @param {number} o.t brace thickness [m]
 * @param {number} o.L chord length [m]
 * @param {number} o.thetaDeg brace angle [deg]
 * @param {number} [o.sigmaX=0] nominal brace axial stress [Pa]
 * @param {number} [o.sigmaMy=0] nominal brace in-plane bending stress [Pa]
 * @param {number} [o.sigmaMz=0] nominal brace out-of-plane bending stress [Pa]
 * @param {'fixed'|'general'} [o.chordEnds='general']
 * @param {number} [o.C=0.7]
 * @returns {{parameters:object, scf:object, chord:object, brace:object, warnings:string[]}} hot-spot results as hotSpotStresses
 */
export function tyJointHotSpots({ D, T, d, t, L, thetaDeg, sigmaX = 0, sigmaMy = 0, sigmaMz = 0, chordEnds = 'general', C = 0.7 }) {
  const p = scfParameters({ D, T, d, t, L, thetaDeg });
  const scf = scfTY({ ...p, chordEnds, C });
  const chord = hotSpotStresses({
    sigmaX, sigmaMy, sigmaMz, scfAC: scf.axial.chordCrown, scfAS: scf.axial.chordSaddle, scfMIP: scf.ipb.chordCrown, scfMOP: scf.opb.chordSaddle,
  });
  const brace = hotSpotStresses({
    sigmaX, sigmaMy, sigmaMz, scfAC: scf.axial.braceCrown, scfAS: scf.axial.braceSaddle, scfMIP: scf.ipb.braceCrown, scfMOP: scf.opb.braceSaddle,
  });
  return { parameters: p, scf, chord, brace, warnings: scf.warnings };
}
