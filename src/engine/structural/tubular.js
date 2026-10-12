/**
 * Strength checks of unstiffened cylindrical tubular members.
 *
 *  - ISO 19902 §13 (partial resistance factor format; equations as in the 2007 edition, Clause 13.2–13.4)
 *  - API RP 2A-LRFD (1st ed., 1993) Section D (without the hydrostatic interaction equations)
 *  - API RP 2A-WSD (21st/22nd ed.) §3.2–3.3 (allowable stress format)
 *
 * Sign convention: axial force N > 0 is tension. Bending moments My, Mz about the two section axes
 * (a tube has no preferred axis; "y" and "z" only pair the moments with the effective length factors).
 * All quantities SI (N, N·m, Pa, m). Every check returns its utilisation with the governing clause.
 *
 * Not covered: conical transitions, grouted members, dented / corroded members, ring-stiffener design,
 * the API RP 2A-LRFD hydrostatic interaction equations.
 */

function pos(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}

function geom(D, t) {
  pos('D', D); pos('t', t);
  if (t > D / 2) throw new RangeError('t must not exceed D/2');
  const d = D - 2 * t;
  const A = Math.PI / 4 * (D * D - d * d);
  const I = Math.PI / 64 * (D ** 4 - d ** 4);
  return { A, I, Ip: 2 * I, Ze: 2 * I / D, Zp: (D ** 3 - d ** 3) / 6, r: Math.sqrt(I / A) };
}

/**
 * Elastic hoop-buckling coefficient Ch as a function of the geometric parameter μ = (L/D)√(2D/t).
 * ISO 19902 §13.2.6.2: Ch = 0.44 t/D (μ ≥ 1.6 D/t); 0.44 t/D + 0.21 (D/t)³/μ⁴ (0.825 D/t ≤ μ < 1.6 D/t);
 * 0.737/(μ − 0.579) (1.5 ≤ μ < 0.825 D/t); 0.80 (μ < 1.5).
 * API RP 2A-WSD §3.2.5.b: as above except 0.736/(M − 0.636) for 3.5 ≤ M < 0.825 D/t and 0.755/(M − 0.559)
 * for 1.5 ≤ M < 3.5.
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.t wall thickness [m]
 * @param {number} o.L length between stiffening rings, diaphragms or end connections [m]
 * @param {'ISO'|'API'} [o.code='ISO'] coefficient set
 * @returns {{Ch:number, mu:number}} coefficient [-] and geometric parameter [-]
 */
export function hoopBucklingCoefficient({ D, t, L, code = 'ISO' }) {
  pos('D', D); pos('t', t); pos('L', L);
  const dt = D / t;
  const mu = (L / D) * Math.sqrt(2 * dt);
  let Ch;
  if (mu >= 1.6 * dt) Ch = 0.44 / dt;
  else if (mu >= 0.825 * dt) Ch = 0.44 / dt + 0.21 * dt ** 3 / mu ** 4;
  else if (code === 'API') {
    if (mu >= 3.5) Ch = 0.736 / (mu - 0.636);
    else if (mu >= 1.5) Ch = 0.755 / (mu - 0.559);
    else Ch = 0.8;
  } else if (mu >= 1.5) Ch = 0.737 / (mu - 0.579);
  else Ch = 0.8;
  return { Ch, mu };
}

/**
 * Representative strengths of a tubular member to ISO 19902 §13.2.
 *  - tension ft = fy (13.2.2)
 *  - local buckling fyc = fy for fy/fxe ≤ 0.170, (1.047 − 0.274 fy/fxe) fy otherwise, fxe = 2 Cx E t/D (13.2.3.3)
 *  - column buckling fc = (1 − 0.278 λ²) fyc for λ ≤ 1.34, 0.9 fyc/λ² otherwise, λ = (K L/(π r))√(fyc/E) (13.2.3.2)
 *  - bending fb = (Zp/Ze) fy for fy D/(E t) ≤ 0.0517; [1.13 − 2.58 fy D/(E t)](Zp/Ze) fy up to 0.1034;
 *    [0.94 − 0.76 fy D/(E t)](Zp/Ze) fy up to 120 fy/E (13.2.4)
 *  - shear fv = fy/√3 (13.2.5)
 *  - hoop buckling fh = fy (fhe > 2.44 fy); 0.7 fy (fhe/fy)^0.4 (0.55 fy < fhe ≤ 2.44 fy); fhe (fhe ≤ 0.55 fy),
 *    fhe = 2 Ch E t/D (13.2.6.2)
 * Validity: D/t ≤ 120, t ≥ 6 mm, fy ≤ 500 MPa (warnings are issued outside).
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.t wall thickness [m]
 * @param {number} o.fy yield strength [Pa]
 * @param {number} [o.E=2.05e11] Young's modulus [Pa]
 * @param {number} [o.L] unbraced length [m] (needed for column buckling)
 * @param {number} [o.K=1] effective length factor [-]
 * @param {number} [o.Lr=o.L] length between rings / diaphragms for hoop buckling [m]
 * @param {number} [o.Cx=0.3] elastic critical buckling coefficient [-]
 * @returns {{A:number, I:number, Ze:number, Zp:number, r:number, ft:number, fxe:number, fyc:number, lambda:number, fc:number, fe:number,
 *   fb:number, fv:number, fhe:number, fh:number, Ch:number, mu:number, Np:number, Mp:number, Pc:number, warnings:string[]}}
 *   areas/moduli (SI), strengths [Pa], λ [-], Np = A fy [N], Mp = Zp fy [N·m], Pc = A fc [N]
 */
export function isoTubularStrengths({ D, t, fy, E = 2.05e11, L, K = 1, Lr, Cx = 0.3 }) {
  const g = geom(D, t);
  pos('fy', fy); pos('E', E);
  const warnings = [];
  if (D / t > 120) warnings.push('D/t exceeds 120: outside the range of ISO 19902 §13.');
  if (t < 0.006) warnings.push('t is below 6 mm: outside the range of ISO 19902 §13.');
  if (fy > 500e6) warnings.push('fy exceeds 500 MPa: outside the range of ISO 19902 §13.');
  const fxe = 2 * Cx * E * t / D;
  const fyc = fy / fxe <= 0.170 ? fy : (1.047 - 0.274 * fy / fxe) * fy;
  let lambda = 0; let fc = fyc; let fe = Infinity;
  if (L !== undefined) {
    pos('L', L); pos('K', K);
    lambda = K * L / (Math.PI * g.r) * Math.sqrt(fyc / E);
    fc = lambda <= 1.34 ? (1 - 0.278 * lambda * lambda) * fyc : 0.9 * fyc / (lambda * lambda);
    fe = Math.PI ** 2 * E / (K * L / g.r) ** 2;
  }
  const x = fy * D / (E * t);
  const shape = g.Zp / g.Ze;
  let fb;
  if (x <= 0.0517) fb = shape * fy;
  else if (x <= 0.1034) fb = (1.13 - 2.58 * x) * shape * fy;
  else {
    fb = (0.94 - 0.76 * x) * shape * fy;
    if (x > 120 * fy / E) warnings.push('fy D/(E t) exceeds 120 fy/E: bending strength formula extrapolated.');
  }
  let fhe = Infinity; let fh = fy; let Ch = NaN; let mu = NaN;
  const lr = Lr === undefined ? L : Lr;
  if (lr !== undefined) {
    ({ Ch, mu } = hoopBucklingCoefficient({ D, t, L: lr, code: 'ISO' }));
    fhe = 2 * Ch * E * t / D;
    if (fhe > 2.44 * fy) fh = fy;
    else if (fhe > 0.55 * fy) fh = 0.7 * fy * (fhe / fy) ** 0.4;
    else fh = fhe;
  }
  return {
    ...g, ft: fy, fxe, fyc, lambda, fc, fe, fb, fv: fy / Math.sqrt(3), fhe, fh, Ch, mu,
    Np: g.A * fy, Mp: g.Zp * fy, Pc: g.A * fc, warnings,
  };
}

function summarise(checks) {
  let gov = null;
  for (const c of checks) if (!gov || c.uc > gov.uc) gov = c;
  return gov;
}

// Largest moment amplification 1/(1 − σc/fe) carried through the beam-column equations. At and beyond the Euler
// stress the member has no second-order equilibrium and the factor is unbounded; it is held at this value so that
// the utilisation stays a large finite number that can be ranked, plotted and iterated on.
const EULER_AMP_CAP = 1000;
const eulerAmp = (ratio) => {
  const d = 1 - ratio;
  return d > 1 / EULER_AMP_CAP ? 1 / d : EULER_AMP_CAP;
};
const EULER_WARNING = 'Axial compression is at or beyond the Euler buckling load: the moment amplification is held at '
  + `${EULER_AMP_CAP} and the beam-column utilisation is indicative only`;

/**
 * ISO 19902 §13 member check for a tubular under axial force, biaxial bending, shear, torsion and (optionally)
 * external hydrostatic pressure. Partial resistance factors default to the ISO values
 * γR,t = 1.05, γR,c = 1.18, γR,b = 1.05, γR,v = 1.05, γR,h = 1.25. Forces are factored design actions.
 *
 * Without pressure (13.3): tension γR,t σt/ft + γR,b σb/fb ≤ 1; compression
 * γR,c σc/fc + (γR,b/fb)√[(Cm,y σb,y/(1 − σc/fe,y))² + (Cm,z σb,z/(1 − σc/fe,z))²] ≤ 1 and γR,c σc/fyc + γR,b σb/fb ≤ 1.
 *
 * With pressure (13.4): σh = p D/(2t), capped-end stress σq = 0.5 σh, B = γR,h σh/fh, η = 5 − 4 fh/fy,
 * ft,h = fy[√(1 + 0.09B² − B^2η) − 0.3B], fb,h = fb[…same…];
 * fc,h = ½ fyc[(1 − 0.278λ²) − 2σq/fyc + √((1 − 0.278λ²)² + 1.12 λ² σq/fyc)] for λ ≤ 1.34/√(1 − 2σq/fyc), else 0.9 fyc/λ²;
 * net tension: γR,t σt,c/ft,h + γR,b σb/fb,h ≤ 1; compression: γR,c σc,c/fyc + γR,b σb/fb,h ≤ 1 and the beam-column
 * equation with fc,h and fb,h; elastic local/hoop interaction (σx − 0.5 fhe/γR,h)/(fxe/γR,c − 0.5 fhe/γR,h) + (γR,h σh/fhe)² ≤ 1
 * when σx > 0.5 fhe/γR,h and fxe > 0.5 fhe.
 *
 * Beyond the Euler load (σc ≥ fe about either axis) the amplification 1/(1 − σc/fe) is held at 1000: every
 * utilisation stays finite (and far above 1 when there is any bending; the axial term alone already exceeds 1),
 * `beyondEuler` is true and a warning is added.
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.t wall thickness [m]
 * @param {number} o.fy yield strength [Pa]
 * @param {number} [o.E=2.05e11] Young's modulus [Pa]
 * @param {number} o.L unbraced length [m]
 * @param {number} [o.Ky=1] effective length factor for bending about y [-]
 * @param {number} [o.Kz=1] effective length factor for bending about z [-]
 * @param {number} [o.Cmy=0.85] moment reduction factor [-]
 * @param {number} [o.Cmz=0.85] moment reduction factor [-]
 * @param {number} [o.N=0] axial force [N], tension positive
 * @param {number} [o.My=0] bending moment [N·m]
 * @param {number} [o.Mz=0] bending moment [N·m]
 * @param {number} [o.V=0] resultant beam shear [N]
 * @param {number} [o.T=0] torsional moment [N·m]
 * @param {number} [o.p=0] factored external hydrostatic pressure [Pa]
 * @param {number} [o.Lr=o.L] ring spacing for hoop buckling [m]
 * @param {boolean} [o.cappedEndIncluded=false] true when N already contains the capped-end (end pressure) force
 * @param {object} [o.gamma] partial resistance factors { t, c, b, v, h }
 * @returns {{stresses:object, strengths:object, checks:Array<{name:string, uc:number, ref:string}>, governing:object,
 *   warnings:string[], beyondEuler:boolean}}
 *   stresses and strengths in Pa, utilisations dimensionless; beyondEuler = compressive stress at
 *   or above the Euler stress about either axis
 */
export function isoTubularCheck({
  D, t, fy, E = 2.05e11, L, Ky = 1, Kz = 1, Cmy = 0.85, Cmz = 0.85, N = 0, My = 0, Mz = 0, V = 0, T = 0, p = 0, Lr,
  cappedEndIncluded = false, gamma = {},
}) {
  const g = { t: 1.05, c: 1.18, b: 1.05, v: 1.05, h: 1.25, ...gamma };
  const Kmax = Math.max(Ky, Kz);
  const S = isoTubularStrengths({ D, t, fy, E, L, K: Kmax, Lr });
  if (!(p >= 0)) throw new RangeError('p must be non-negative');
  const sa = N / S.A;
  const sby = Math.abs(My) / S.Ze;
  const sbz = Math.abs(Mz) / S.Ze;
  const sb = Math.hypot(sby, sbz);
  const tauB = 2 * Math.abs(V) / S.A;
  const tauT = Math.abs(T) * D / (2 * S.Ip);
  const sh = p * D / (2 * t);
  const sq = 0.5 * sh;
  const fey = Math.PI ** 2 * E / (Ky * L / S.r) ** 2;
  const fez = Math.PI ** 2 * E / (Kz * L / S.r) ** 2;
  const checks = [];
  const warnings = [...S.warnings];
  const add = (name, uc, ref) => checks.push({ name, uc, ref });
  // axial stress without capped-end action (tension positive)
  const saEx = cappedEndIncluded ? sa + sq : sa;
  if (saEx >= 0) add('Axial tension', g.t * saEx / S.ft, 'ISO 19902 §13.2.2');
  else add('Axial compression (column buckling)', g.c * -saEx / S.fc, 'ISO 19902 §13.2.3');
  add('Bending', g.b * sb / S.fb, 'ISO 19902 §13.2.4');
  add('Beam shear', g.v * tauB / S.fv, 'ISO 19902 §13.2.5.1');
  add('Torsional shear', g.v * tauT / S.fv, 'ISO 19902 §13.2.5.2');
  const amp = (sc, fe) => eulerAmp(sc / fe);
  const beyondEuler = saEx < 0 && -saEx >= Math.min(fey, fez);
  if (beyondEuler) warnings.push(EULER_WARNING);
  if (sh === 0) {
    if (saEx >= 0) {
      add('Tension + bending', g.t * saEx / S.ft + g.b * sb / S.fb, 'ISO 19902 §13.3.2');
    } else {
      const sc = -saEx;
      add('Compression + bending (beam-column)',
        g.c * sc / S.fc + g.b / S.fb * Math.hypot(Cmy * sby * amp(sc, fey), Cmz * sbz * amp(sc, fez)), 'ISO 19902 §13.3.3 (stability)');
      add('Compression + bending (local strength)', g.c * sc / S.fyc + g.b * sb / S.fb, 'ISO 19902 §13.3.3 (strength)');
    }
  } else {
    add('Hoop buckling', g.h * sh / S.fh, 'ISO 19902 §13.2.6');
    const B = Math.min(g.h * sh / S.fh, 1);
    const eta = 5 - 4 * S.fh / fy;
    const red = Math.max(Math.sqrt(Math.max(1 + 0.09 * B * B - B ** (2 * eta), 0)) - 0.3 * B, 0);
    const fth = fy * red;
    const fbh = S.fb * red;
    const net = saEx - sq; // axial stress including the capped-end action
    const ratio = (num, den) => (num === 0 ? 0 : (den > 0 ? num / den : Infinity));
    let scc = 0;
    if (net >= 0) {
      add('Tension + bending + hydrostatic pressure', ratio(g.t * net, fth) + ratio(g.b * sb, fbh), 'ISO 19902 §13.4.2');
    } else {
      scc = -net;
      add('Compression + bending + hydrostatic pressure (local strength)', g.c * scc / S.fyc + ratio(g.b * sb, fbh), 'ISO 19902 §13.4.3');
    }
    if (saEx < 0) {
      const sc = -saEx;
      const lam = S.lambda;
      const a = sq / S.fyc;
      let fch;
      if (a < 0.5 && lam <= 1.34 / Math.sqrt(1 - 2 * a)) {
        const c0 = 1 - 0.278 * lam * lam;
        fch = 0.5 * S.fyc * (c0 - 2 * a + Math.sqrt(c0 * c0 + 1.12 * lam * lam * a));
      } else {
        fch = 0.9 * S.fyc / (lam * lam);
      }
      add('Compression + bending + hydrostatic pressure (beam-column)',
        ratio(g.c * sc, fch) + ratio(g.b * Math.hypot(Cmy * sby * amp(sc, fey), Cmz * sbz * amp(sc, fez)), fbh), 'ISO 19902 §13.4.3 (stability)');
    }
    const sx = scc + sb;
    if (sx > 0.5 * S.fhe / g.h && S.fxe > 0.5 * S.fhe) {
      add('Elastic local / hoop buckling interaction',
        (sx - 0.5 * S.fhe / g.h) / (S.fxe / g.c - 0.5 * S.fhe / g.h) + (g.h * sh / S.fhe) ** 2, 'ISO 19902 §13.4.3 (elastic interaction)');
    }
  }
  return {
    stresses: { axial: sa, axialExcludingCappedEnd: saEx, bendingY: sby, bendingZ: sbz, bending: sb, shear: tauB, torsion: tauT, hoop: sh, cappedEnd: sq },
    strengths: S, checks, governing: summarise(checks), warnings, beyondEuler,
  };
}

/**
 * API RP 2A-LRFD (1993) Section D member check without hydrostatic interaction. Resistance factors
 * φt = 0.95, φc = 0.85, φb = 0.95, φv = 0.95, φh = 0.80.
 *  - D.2.2 column: Fcn = (1 − 0.25 λ²) Fy for λ < √2, Fy/λ² otherwise, λ = (K L/(π r))√(Fy/E); for D/t > 60, Fy is
 *    replaced by the local buckling stress min(Fxe, Fxc), Fxe = 2 Cx E t/D (Cx = 0.3), Fxc = [1.64 − 0.23 (D/t)^¼] Fy
 *  - D.2.3 bending: as ISO 19902 §13.2.4 (same three-range expression in Fy D/(E t))
 *  - D.2.4 shear Fvn = Fy/√3 with fv = 2V/A; torsion fvt = T D/(2 Ip)
 *  - D.2.5 hoop: Fhc = Fhe (Fhe ≤ 0.55 Fy), 0.7 Fy (Fhe/Fy)^0.4 ≤ Fy otherwise; fh ≤ φh Fhc
 *  - D.3.1 tension + bending: 1 − cos(π ft/(2 φt Fy)) + √(fby² + fbz²)/(φb Fbn) ≤ 1
 *  - D.3.2 compression + bending: fc/(φc Fcn) + (1/(φb Fbn))√[(Cmy fby/(1 − fc/(φc Fey)))² + (Cmz fbz/(1 − fc/(φc Fez)))²] ≤ 1
 *    and 1 − cos(π fc/(2 φc Fxc)) + √(fby² + fbz²)/(φb Fbn) ≤ 1
 * @param {object} o same geometry / force arguments as isoTubularCheck (p only used for the hoop check)
 * @param {object} [o.phi] resistance factors { t, c, b, v, h }
 * @returns {{stresses:object, strengths:object, checks:Array, governing:object, warnings:string[], beyondEuler:boolean}}
 *   stresses and strengths in Pa; beyondEuler = compression at or above the (factored) Euler stress, where the
 *   moment amplification is held at 1000 so that the beam-column utilisation stays finite
 */
export function apiLrfdTubularCheck({
  D, t, fy, E = 2.0e11, L, Ky = 1, Kz = 1, Cmy = 0.85, Cmz = 0.85, N = 0, My = 0, Mz = 0, V = 0, T = 0, p = 0, Lr, phi = {},
}) {
  const f = { t: 0.95, c: 0.85, b: 0.95, v: 0.95, h: 0.80, ...phi };
  const g = geom(D, t);
  pos('fy', fy); pos('E', E); pos('L', L);
  const warnings = [];
  let beyondEuler = false;
  const dt = D / t;
  if (dt > 300) warnings.push('D/t exceeds 300: outside the range of API RP 2A-LRFD Section D.');
  const Fxe = 2 * 0.3 * E * t / D;
  const Fxc = dt <= 60 ? fy : Math.min((1.64 - 0.23 * dt ** 0.25) * fy, fy);
  const Fyl = dt <= 60 ? fy : Math.min(Fxe, Fxc);
  const K = Math.max(Ky, Kz);
  const lambda = K * L / (Math.PI * g.r) * Math.sqrt(Fyl / E);
  const Fcn = lambda < Math.SQRT2 ? (1 - 0.25 * lambda * lambda) * Fyl : Fyl / (lambda * lambda);
  const x = fy * D / (E * t);
  const shape = g.Zp / g.Ze;
  const Fbn = x <= 0.0517 ? shape * fy : (x <= 0.1034 ? (1.13 - 2.58 * x) * shape * fy : (0.94 - 0.76 * x) * shape * fy);
  const Fvn = fy / Math.sqrt(3);
  const Fey = Math.PI ** 2 * E / (Ky * L / g.r) ** 2;
  const Fez = Math.PI ** 2 * E / (Kz * L / g.r) ** 2;
  const sa = N / g.A;
  const sby = Math.abs(My) / g.Ze; const sbz = Math.abs(Mz) / g.Ze; const sb = Math.hypot(sby, sbz);
  const checks = [];
  const add = (name, uc, ref) => checks.push({ name, uc, ref });
  add('Bending', sb / (f.b * Fbn), 'API RP 2A-LRFD D.2.3');
  add('Beam shear', 2 * Math.abs(V) / g.A / (f.v * Fvn), 'API RP 2A-LRFD D.2.4.1');
  add('Torsional shear', Math.abs(T) * D / (2 * g.Ip) / (f.v * Fvn), 'API RP 2A-LRFD D.2.4.2');
  if (sa >= 0) {
    add('Axial tension', sa / (f.t * fy), 'API RP 2A-LRFD D.2.1');
    add('Tension + bending', 1 - Math.cos(Math.PI / 2 * Math.min(sa / (f.t * fy), 2)) + sb / (f.b * Fbn), 'API RP 2A-LRFD D.3.1');
  } else {
    const fc = -sa;
    add('Axial compression', fc / (f.c * Fcn), 'API RP 2A-LRFD D.2.2');
    const amp = (fe) => eulerAmp(fc / (f.c * fe));
    if (fc >= f.c * Math.min(Fey, Fez)) { beyondEuler = true; warnings.push(EULER_WARNING); }
    add('Compression + bending (beam-column)', fc / (f.c * Fcn) + Math.hypot(Cmy * sby * amp(Fey), Cmz * sbz * amp(Fez)) / (f.b * Fbn), 'API RP 2A-LRFD D.3.2-1');
    add('Compression + bending (local strength)', 1 - Math.cos(Math.PI / 2 * Math.min(fc / (f.c * Fxc), 2)) + sb / (f.b * Fbn), 'API RP 2A-LRFD D.3.2-2');
  }
  let Fhe = Infinity; let Fhc = fy; let sh = 0;
  if (p > 0) {
    const { Ch } = hoopBucklingCoefficient({ D, t, L: Lr === undefined ? L : Lr, code: 'API' });
    Fhe = 2 * Ch * E * t / D;
    Fhc = Fhe <= 0.55 * fy ? Fhe : Math.min(0.7 * fy * (Fhe / fy) ** 0.4, fy);
    sh = p * D / (2 * t);
    add('Hoop buckling', sh / (f.h * Fhc), 'API RP 2A-LRFD D.2.5');
    warnings.push('Hydrostatic interaction with axial force and bending (API RP 2A-LRFD D.3.3/D.3.4) is not evaluated.');
  }
  return {
    stresses: { axial: sa, bendingY: sby, bendingZ: sbz, bending: sb, hoop: sh },
    strengths: { ...g, Fxe, Fxc, lambda, Fcn, Fbn, Fvn, Fey, Fez, Fhe, Fhc }, checks, governing: summarise(checks), warnings, beyondEuler,
  };
}

/**
 * API RP 2A-WSD (21st/22nd ed.) §3.2–3.3 allowable-stress check of a tubular member.
 *  - 3.2.1 tension Ft = 0.6 Fy
 *  - 3.2.2 compression (AISC): Fa = [1 − (Kl/r)²/(2Cc²)]Fy / [5/3 + 3(Kl/r)/(8Cc) − (Kl/r)³/(8Cc³)] for Kl/r < Cc,
 *    12π²E/(23 (Kl/r)²) otherwise, Cc = √(2π²E/Fy); for D/t > 60, Fy is replaced by min(Fxe, Fxc) with
 *    Fxe = 2 C E t/D (C = 0.3), Fxc = Fy[1.64 − 0.23 (D/t)^¼] ≤ Fxe
 *  - 3.2.3 bending Fb = 0.75 Fy (D/t ≤ 10340/Fy); [0.84 − 1.74 Fy D/(E t)]Fy (≤ 20680/Fy); [0.72 − 0.58 Fy D/(E t)]Fy (≤ 300), Fy in MPa
 *  - 3.2.4 shear fv = V/(0.5A) ≤ 0.4 Fy; torsion fvt = T D/(2 Ip) ≤ 0.4 Fy
 *  - 3.2.5 hoop: fh = p D/(2t) ≤ Fhc/SFh; Fhc = Fhe (≤ 0.55Fy), 0.45Fy + 0.18Fhe (≤ 1.6Fy), 1.31Fy/(1.15 + Fy/Fhe) (≤ 6.2Fy), Fy
 *  - 3.3.1 compression + bending: Eq. 3.3.1-1 (per-axis amplification form 3.3.1-4), 3.3.1-2, and 3.3.1-3 when fa/Fa ≤ 0.15;
 *    tension + bending: fa/(0.6Fy) + fb/Fb ≤ 1
 *  - 3.3.3 tension + hydrostatic: A² + B² + 2ν|A|B ≤ 1, A = (fa + fb − 0.5fh) SFx/Fy, B = fh SFh/Fhc, ν = 0.3
 *  - 3.3.4 compression + hydrostatic: (fa + 0.5fh) SFx/Fxc + fb SFb/Fy ≤ 1 (3.3.4-1); (fx − 0.5Fha)/(Faa − 0.5Fha) + (fh/Fha)² ≤ 1
 *    for fx > 0.5 Fha (3.3.4-3), Faa = Fxe/SFx, Fha = Fhe/SFh, fx = fa + fb + 0.5fh
 * The one-third allowable stress increase for extreme conditions is applied through `allowableFactor` (1.0 or 1.333).
 * @param {object} o geometry / forces as isoTubularCheck (unfactored loads); N excludes the capped-end force
 * @param {number} [o.allowableFactor=1] multiplier on basic allowable stresses [-]
 * @param {number} [o.SFxTension=1.67] safety factor, axial tension with hydrostatic pressure [-]
 * @param {number} [o.SFxCompression=2.0] safety factor, axial compression with hydrostatic pressure [-] (1.67–2.0)
 * @param {number} [o.SFh=2.0] safety factor, hoop compression [-]
 * @returns {{stresses:object, strengths:object, checks:Array, governing:object, warnings:string[], beyondEuler:boolean}}
 *   stresses and strengths in Pa; beyondEuler = compression at or above the (factored) Euler stress, where the
 *   moment amplification is held at 1000 so that the beam-column utilisation stays finite
 */
export function apiWsdTubularCheck({
  D, t, fy, E = 2.0e11, L, Ky = 1, Kz = 1, Cmy = 0.85, Cmz = 0.85, N = 0, My = 0, Mz = 0, V = 0, T = 0, p = 0, Lr,
  allowableFactor = 1, SFxTension = 1.67, SFxCompression = 2.0, SFh = 2.0,
}) {
  const g = geom(D, t);
  pos('fy', fy); pos('E', E); pos('L', L); pos('allowableFactor', allowableFactor);
  const warnings = [];
  let beyondEuler = false;
  const dt = D / t;
  if (dt > 300) warnings.push('D/t exceeds 300: outside the range of API RP 2A-WSD §3.2.');
  const af = allowableFactor;
  const Fxe = 2 * 0.3 * E * t / D;
  const Fxc = dt <= 60 ? fy : Math.min(fy * (1.64 - 0.23 * dt ** 0.25), Fxe);
  const Fyl = dt <= 60 ? fy : Math.min(Fxe, Fxc);
  const Fa = (K) => {
    const kl = K * L / g.r;
    const Cc = Math.sqrt(2 * Math.PI ** 2 * E / Fyl);
    if (kl < Cc) {
      return (1 - kl * kl / (2 * Cc * Cc)) * Fyl / (5 / 3 + 3 * kl / (8 * Cc) - kl ** 3 / (8 * Cc ** 3));
    }
    return 12 * Math.PI ** 2 * E / (23 * kl * kl);
  };
  const FaMin = af * Fa(Math.max(Ky, Kz));
  const fyMPa = fy / 1e6;
  const x = fy * D / (E * t);
  let Fb;
  if (dt <= 10340 / fyMPa) Fb = 0.75 * fy;
  else if (dt <= 20680 / fyMPa) Fb = (0.84 - 1.74 * x) * fy;
  else Fb = (0.72 - 0.58 * x) * fy;
  const FbBasic = Fb;
  Fb *= af;
  const Ft = 0.6 * fy * af;
  const Fv = 0.4 * fy * af;
  const Fe = (K) => af * 12 * Math.PI ** 2 * E / (23 * (K * L / g.r) ** 2);
  const fa = N / g.A;
  const fby = Math.abs(My) / g.Ze; const fbz = Math.abs(Mz) / g.Ze; const fb = Math.hypot(fby, fbz);
  const checks = [];
  const add = (name, uc, ref) => checks.push({ name, uc, ref });
  add('Bending', fb / Fb, 'API RP 2A-WSD 3.2.3');
  add('Beam shear', Math.abs(V) / (0.5 * g.A) / Fv, 'API RP 2A-WSD 3.2.4a');
  add('Torsional shear', Math.abs(T) * D / (2 * g.Ip) / Fv, 'API RP 2A-WSD 3.2.4b');
  if (fa >= 0) {
    add('Axial tension', fa / Ft, 'API RP 2A-WSD 3.2.1');
    add('Tension + bending', fa / Ft + fb / Fb, 'API RP 2A-WSD 3.3.1b');
  } else {
    const fc = -fa;
    add('Axial compression', fc / FaMin, 'API RP 2A-WSD 3.2.2');
    add('Compression + bending (strength)', fc / (0.6 * Fyl * af) + fb / Fb, 'API RP 2A-WSD Eq. 3.3.1-2');
    if (fc / FaMin <= 0.15) {
      add('Compression + bending (low axial)', fc / FaMin + fb / Fb, 'API RP 2A-WSD Eq. 3.3.1-3');
    } else {
      const amp = (K) => eulerAmp(fc / Fe(K));
      if (fc >= Math.min(Fe(Ky), Fe(Kz))) { beyondEuler = true; warnings.push(EULER_WARNING); }
      add('Compression + bending (stability)', fc / FaMin + Math.hypot(Cmy * fby * amp(Ky), Cmz * fbz * amp(Kz)) / Fb, 'API RP 2A-WSD Eq. 3.3.1-4');
    }
  }
  let Fhe = Infinity; let Fhc = fy; let fh = 0;
  if (p > 0) {
    const { Ch } = hoopBucklingCoefficient({ D, t, L: Lr === undefined ? L : Lr, code: 'API' });
    Fhe = 2 * Ch * E * t / D;
    if (Fhe <= 0.55 * fy) Fhc = Fhe;
    else if (Fhe <= 1.6 * fy) Fhc = 0.45 * fy + 0.18 * Fhe;
    else if (Fhe <= 6.2 * fy) Fhc = 1.31 * fy / (1.15 + fy / Fhe);
    else Fhc = fy;
    fh = p * D / (2 * t);
    add('Hoop buckling', fh * SFh / Fhc, 'API RP 2A-WSD Eq. 3.2.5-1');
    const SFb = fy / FbBasic;
    if (fa + 0 >= 0 && fa + fb - 0.5 * fh >= 0) {
      const A = (fa + fb - 0.5 * fh) / fy * SFxTension;
      const B = fh / Fhc * SFh;
      add('Tension + hydrostatic pressure', A * A + B * B + 2 * 0.3 * Math.abs(A) * B, 'API RP 2A-WSD Eq. 3.3.3-1');
    } else {
      const fc = Math.max(-fa, 0);
      add('Compression + hydrostatic pressure', (fc + 0.5 * fh) / Fxc * SFxCompression + fb / fy * SFb, 'API RP 2A-WSD Eq. 3.3.4-1');
      const Faa = Fxe / SFxCompression;
      const Fha = Fhe / SFh;
      const fx = fc + fb + 0.5 * fh;
      if (fx > 0.5 * Fha && Faa > 0.5 * Fha) {
        add('Elastic local / hoop buckling interaction', (fx - 0.5 * Fha) / (Faa - 0.5 * Fha) + (fh / Fha) ** 2, 'API RP 2A-WSD Eq. 3.3.4-3');
      }
    }
  }
  return {
    stresses: { axial: fa, bendingY: fby, bendingZ: fbz, bending: fb, hoop: fh },
    strengths: { ...g, Ft, Fa: FaMin, Fb, Fv, Fxe, Fxc, Fhe, Fhc }, checks, governing: summarise(checks), warnings, beyondEuler,
  };
}

/**
 * Code-check a tubular frame member from its section forces at one station.
 * @param {object} o
 * @param {'ISO'|'API-LRFD'|'API-WSD'} [o.code='ISO'] design code
 * @param {object} o.section tubeSection record or { D, t } [m]
 * @param {number} o.fy yield strength [Pa]
 * @param {object} o.forces section forces { N, Vy, Vz, T, My, Mz } [N, N·m]
 * @param {number} o.L unbraced length [m]
 * @param {object} [o.options] further arguments of the selected check (Ky, Kz, Cmy, Cmz, E, p, Lr, gamma, …)
 * @returns {object} result of isoTubularCheck / apiLrfdTubularCheck / apiWsdTubularCheck
 */
export function tubularMemberCheck({ code = 'ISO', section, fy, forces, L, options = {} }) {
  if (!section || !(section.D > 0)) throw new RangeError('section must be a tube with D and t');
  const f = forces || {};
  const args = {
    ...options, D: section.D, t: section.t, fy, L, N: f.N || 0, My: f.My || 0, Mz: f.Mz || 0,
    V: Math.hypot(f.Vy || 0, f.Vz || 0), T: f.T || 0,
  };
  if (code === 'ISO') return isoTubularCheck(args);
  if (code === 'API-LRFD') return apiLrfdTubularCheck(args);
  if (code === 'API-WSD') return apiWsdTubularCheck(args);
  throw new RangeError(`code "${code}" is not recognised`);
}
