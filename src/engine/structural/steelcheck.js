/**
 * Member and connection checks for open (and hollow) steel sections.
 *
 *  - EN 1993-1-1:2005 — cross-section classification (Table 5.2), resistances (§6.2), flexural buckling (§6.3.1),
 *    lateral-torsional buckling general case (§6.3.2.2), beam-column interaction Method 2 (Annex B, Tables B.1–B.3)
 *    and the conservative linear cross-section summation of §6.2.1(7).
 *  - ANSI/AISC 360-16 — E3 flexural buckling, F2 major-axis flexure of doubly symmetric compact I-shapes,
 *    F6 minor-axis flexure, H1 interaction (LRFD or ASD).
 *  - EN 1993-1-8:2005 — bolt resistances (Table 3.4), simplified fillet weld method (§4.5.3.3), elastic bolt-group and
 *    weld-group analysis.
 *  - Lifting attachments: padeye (bearing, tear-out, net section, cheek plates) and trunnion stress checks on a
 *    first-principles basis with allowable-stress ratios supplied by the caller.
 *
 * SI units (m, N, N·m, Pa). Axial force N > 0 is tension. Section axes: y strong, z weak (see sections.js).
 * Not covered: class 4 effective sections, torsional / flexural-torsional buckling, plate-girder shear buckling,
 * EN 1993-1-1 Annex A (Method 1), AISC slender-element reductions.
 */

function pos(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}

/* ------------------------------------------------------------------ EN 1993-1-1 */

/**
 * Imperfection factor α of the EN 1993-1-1 buckling curves (Table 6.1): a0 0.13, a 0.21, b 0.34, c 0.49, d 0.76.
 * @param {'a0'|'a'|'b'|'c'|'d'} curve
 * @returns {number} α [-]
 */
export function ec3ImperfectionFactor(curve) {
  const a = { a0: 0.13, a: 0.21, b: 0.34, c: 0.49, d: 0.76 }[curve];
  if (a === undefined) throw new RangeError(`curve "${curve}" is not one of a0, a, b, c, d`);
  return a;
}

/**
 * Reduction factor χ for flexural buckling (EN 1993-1-1 Eq. 6.49): χ = 1/(Φ + √(Φ² − λ̄²)) ≤ 1,
 * Φ = 0.5[1 + α(λ̄ − 0.2) + λ̄²]; χ = 1 for λ̄ ≤ 0.2.
 * @param {number} lambdaBar non-dimensional slenderness λ̄ [-]
 * @param {'a0'|'a'|'b'|'c'|'d'} curve buckling curve
 * @returns {number} χ [-]
 */
export function ec3Chi(lambdaBar, curve) {
  if (!(lambdaBar >= 0)) throw new RangeError('lambdaBar must be non-negative');
  const alpha = ec3ImperfectionFactor(curve);
  if (lambdaBar <= 0.2) return 1;
  const phi = 0.5 * (1 + alpha * (lambdaBar - 0.2) + lambdaBar * lambdaBar);
  return Math.min(1, 1 / (phi + Math.sqrt(phi * phi - lambdaBar * lambdaBar)));
}

/**
 * Buckling curves for an I/H section per EN 1993-1-1 Table 6.2 (flexural) and Table 6.4 (lateral-torsional,
 * general case). Rolled: h/b > 1.2 → y: a, z: b (tf ≤ 40 mm), y: b, z: c (40 < tf ≤ 100 mm); h/b ≤ 1.2 → y: b, z: c
 * (tf ≤ 100 mm), d, d above; S460: a0 instead of a/b for tf ≤ 40 mm (h/b > 1.2: a0, a0; h/b ≤ 1.2 and tf ≤ 100: a, a),
 * c, c above 100 mm. Welded: y: b, z: c (tf ≤ 40 mm), y: c, z: d above. LTB: rolled h/b ≤ 2 → a, > 2 → b; welded
 * h/b ≤ 2 → c, > 2 → d.
 * @param {object} o
 * @param {number} o.h depth [m]
 * @param {number} o.b width [m]
 * @param {number} o.tf flange thickness [m]
 * @param {boolean} [o.welded=false]
 * @param {number} [o.fy=355e6] yield strength [Pa] (≥ 460 MPa selects the S460 column)
 * @returns {{y:string, z:string, lt:string}} curve names
 */
export function ec3BucklingCurvesI({ h, b, tf, welded = false, fy = 355e6 }) {
  pos('h', h); pos('b', b); pos('tf', tf);
  const s460 = fy >= 460e6;
  const hb = h / b;
  let y; let z;
  if (welded) {
    if (tf <= 0.040) { y = 'b'; z = 'c'; } else { y = 'c'; z = 'd'; }
  } else if (hb > 1.2) {
    if (tf <= 0.040) { y = s460 ? 'a0' : 'a'; z = s460 ? 'a0' : 'b'; }
    else if (tf <= 0.100) { y = s460 ? 'a' : 'b'; z = s460 ? 'a' : 'c'; }
    else { y = s460 ? 'c' : 'd'; z = s460 ? 'c' : 'd'; }
  } else if (tf <= 0.100) { y = s460 ? 'a' : 'b'; z = s460 ? 'a' : 'c'; }
  else { y = s460 ? 'c' : 'd'; z = s460 ? 'c' : 'd'; }
  const lt = welded ? (hb <= 2 ? 'c' : 'd') : (hb <= 2 ? 'a' : 'b');
  return { y, z, lt };
}

/**
 * Cross-section class of an I/H section or a circular hollow section (EN 1993-1-1 Table 5.2), ε = √(235/fy).
 * Outstand flange in compression: c/t ≤ 9ε, 10ε, 14ε. Web (internal part): pure bending 72ε, 83ε, 124ε; pure compression
 * 33ε, 38ε, 42ε; bending + compression with plastic stress-block parameter α: class 1 396ε/(13α − 1) (α > 0.5) or
 * 36ε/α; class 2 456ε/(13α − 1) or 41.5ε/α; class 3 with elastic stress ratio ψ: 42ε/(0.67 + 0.33ψ) (ψ > −1) or
 * 62ε(1 − ψ)√(−ψ). Circular hollow sections: d/t ≤ 50ε², 70ε², 90ε².
 * For combined loading the web is classified for the given design axial force NEd (compression positive here).
 * @param {object} o
 * @param {object} o.section iSection / weldedISection / tubeSection record
 * @param {number} o.fy yield strength [Pa]
 * @param {number} [o.NEd=0] design axial compression [N] (0 → web in pure bending; Infinity → pure compression)
 * @returns {{class:number, flangeClass:number, webClass:number, epsilon:number, flangeRatio:number, webRatio:number}} classes 1–4 and c/t ratios
 */
export function ec3Classify({ section, fy, NEd = 0 }) {
  pos('fy', fy);
  const eps = Math.sqrt(235e6 / fy);
  const s = section;
  if (!s) throw new RangeError('section is required');
  const classOf = (ratio, limits) => (ratio <= limits[0] ? 1 : (ratio <= limits[1] ? 2 : (ratio <= limits[2] ? 3 : 4)));
  if (s.type === 'tube') {
    const r = s.D / s.t;
    const c = classOf(r, [50 * eps * eps, 70 * eps * eps, 90 * eps * eps]);
    return { class: c, flangeClass: c, webClass: c, epsilon: eps, flangeRatio: r, webRatio: r };
  }
  if (!(s.tw > 0 && s.tf > 0 && s.b > 0 && s.h > 0)) throw new RangeError('section must be an I-section or a tube');
  const r = s.r || 0;
  const cf = (s.b - s.tw - 2 * r) / 2;
  const cw = s.h - 2 * s.tf - 2 * r;
  const fr = cf / s.tf;
  const wr = cw / s.tw;
  const flangeClass = classOf(fr, [9 * eps, 10 * eps, 14 * eps]);
  let webLimits;
  if (NEd === Infinity) webLimits = [33 * eps, 38 * eps, 42 * eps];
  else if (!(NEd > 0)) webLimits = [72 * eps, 83 * eps, 124 * eps];
  else {
    // plastic distribution: α = ½(1 + NEd/(fy tw c)) ≤ 1; elastic: ψ = 2 NEd/(A fy) − 1
    const alpha = Math.min(1, 0.5 * (1 + NEd / (fy * s.tw * cw)));
    const psi = Math.min(1, 2 * NEd / (s.A * fy) - 1);
    const l1 = alpha > 0.5 ? 396 * eps / (13 * alpha - 1) : 36 * eps / alpha;
    const l2 = alpha > 0.5 ? 456 * eps / (13 * alpha - 1) : 41.5 * eps / alpha;
    const l3 = psi > -1 ? 42 * eps / (0.67 + 0.33 * psi) : 62 * eps * (1 - psi) * Math.sqrt(-psi);
    webLimits = [l1, l2, l3];
  }
  const webClass = classOf(wr, webLimits);
  return { class: Math.max(flangeClass, webClass), flangeClass, webClass, epsilon: eps, flangeRatio: fr, webRatio: wr };
}

/**
 * Elastic critical moment for lateral-torsional buckling of a doubly symmetric section loaded at the shear
 * centre: Mcr = C1 (π² E Iz/(k L)²) √[(k/kw)² Iw/Iz + (k L)² G It/(π² E Iz)] (ENV 1993-1-1 Annex F / NCCI SN003, with
 * zg = 0).
 * @param {object} o
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} o.G shear modulus [Pa]
 * @param {number} o.Iz weak-axis second moment [m⁴]
 * @param {number} o.It torsion constant [m⁴]
 * @param {number} o.Iw warping constant [m⁶]
 * @param {number} o.L length between lateral restraints [m]
 * @param {number} [o.C1=1] moment-gradient factor [-] (1.0 uniform moment, 1.127 uniform load, 1.348 mid-span point load — pinned ends)
 * @param {number} [o.k=1] effective length factor for end rotation on plan [-]
 * @param {number} [o.kw=1] effective length factor for end warping [-]
 * @returns {number} Mcr [N·m]
 */
export function elasticCriticalMoment({ E, G, Iz, It, Iw, L, C1 = 1, k = 1, kw = 1 }) {
  pos('E', E); pos('G', G); pos('Iz', Iz); pos('L', L); pos('C1', C1);
  if (!(It >= 0) || !(Iw >= 0)) throw new RangeError('It and Iw must be non-negative');
  const kl = k * L;
  const pe = Math.PI ** 2 * E * Iz / (kl * kl);
  return C1 * pe * Math.sqrt((k / kw) ** 2 * Iw / Iz + G * It / pe);
}

/**
 * Equivalent uniform moment factor for a linear moment diagram, EN 1993-1-1 Table B.3: Cm = 0.6 + 0.4ψ ≥ 0.4.
 * @param {number} psi end-moment ratio −1 ≤ ψ ≤ 1 [-]
 * @returns {number} Cm [-]
 */
export function ec3CmLinear(psi) {
  if (!(psi >= -1 && psi <= 1)) throw new RangeError('psi must lie in [-1, 1]');
  return Math.max(0.4, 0.6 + 0.4 * psi);
}

/**
 * EN 1993-1-1 member check of a doubly symmetric I/H section (class 1–3) under axial force, biaxial bending and
 * major-axis shear.
 *  - §6.2.3/6.2.4 Npl,Rd = A fy/γM0; §6.2.5 Mc,Rd = W fy/γM0 (Wpl for class 1–2, Wel for class 3); §6.2.6 Vpl,Rd = Av fy/(√3 γM0)
 *  - §6.2.8 high shear: for VEd > 0.5 Vpl,Rd the major-axis bending resistance is reduced with ρ = (2VEd/Vpl,Rd − 1)² acting on
 *    the web (My,V,Rd = [Wpl,y − ρ Aw²/(4 tw)] fy/γM0, Aw = hw tw)
 *  - §6.2.1(7) conservative cross-section sum NEd/NRd + My,Ed/My,Rd + Mz,Ed/Mz,Rd ≤ 1
 *  - §6.3.1 flexural buckling Nb,Rd = χ A fy/γM1 about both axes (curves from Table 6.2)
 *  - §6.3.2.2 lateral-torsional buckling, general case: χLT from λ̄LT = √(Wy fy/Mcr) with αLT of Table 6.3/6.4
 *  - §6.3.3 Eq. (6.61), (6.62) with interaction factors of Annex B (Method 2): Table B.2 (members susceptible to torsional
 *    deformations) by default, Table B.1 when `torsionallyRestrained` is true.
 * Tension members are checked for cross-section resistance and (with bending) LTB only.
 * @param {object} o
 * @param {object} o.section iSection record (A, Iy, Iz, J, Iw, Wely, Welz, Wply, Wplz, Asz, h, b, tw, tf, r)
 * @param {number} o.fy yield strength [Pa]
 * @param {number} [o.E=2.1e11] [Pa]
 * @param {number} [o.G=8.1e10] [Pa]
 * @param {number} [o.NEd=0] axial force [N], tension positive
 * @param {number} [o.MyEd=0] maximum major-axis moment [N·m]
 * @param {number} [o.MzEd=0] maximum minor-axis moment [N·m]
 * @param {number} [o.VzEd=0] shear along the web [N]
 * @param {number} o.Lcry buckling length about y [m]
 * @param {number} [o.Lcrz=o.Lcry] buckling length about z [m]
 * @param {number} [o.LLT=o.Lcrz] length between lateral restraints [m]
 * @param {number} [o.C1=1] moment-gradient factor for Mcr [-]
 * @param {number} [o.Cmy=1] [-]
 * @param {number} [o.Cmz=1] [-]
 * @param {number} [o.CmLT=1] [-]
 * @param {boolean} [o.torsionallyRestrained=false] use Table B.1 (no LTB)
 * @param {number} [o.gammaM0=1.0] [-]
 * @param {number} [o.gammaM1=1.0] [-]
 * @returns {{classification:object, curves:object, resistances:object, factors:object, checks:Array<{name:string, uc:number, ref:string}>,
 *   governing:object, warnings:string[]}} resistances [N, N·m], utilisations [-]
 */
export function ec3MemberCheck({
  section, fy, E = 2.1e11, G = 8.1e10, NEd = 0, MyEd = 0, MzEd = 0, VzEd = 0, Lcry, Lcrz, LLT, C1 = 1, Cmy = 1, Cmz = 1, CmLT = 1,
  torsionallyRestrained = false, gammaM0 = 1.0, gammaM1 = 1.0,
}) {
  const s = section;
  pos('fy', fy); pos('Lcry', Lcry);
  if (!s || !(s.A > 0) || !(s.tw > 0)) throw new RangeError('section must be an I-section record');
  const lz = Lcrz === undefined ? Lcry : pos('Lcrz', Lcrz);
  const llt = LLT === undefined ? lz : pos('LLT', LLT);
  const warnings = [];
  const comp = Math.max(-NEd, 0);
  const cls = ec3Classify({ section: s, fy, NEd: comp });
  if (cls.class === 4) warnings.push('Class 4 section: effective properties are not evaluated; elastic gross properties are used.');
  const plastic = cls.class <= 2;
  const Wy = plastic ? s.Wply : s.Wely;
  const Wz = plastic ? s.Wplz : s.Welz;
  const NRk = s.A * fy; const MyRk = Wy * fy; const MzRk = Wz * fy;
  const NplRd = NRk / gammaM0;
  const Av = s.Asz;
  const VplRd = Av * fy / (Math.sqrt(3) * gammaM0);
  let MyRd = MyRk / gammaM0;
  const MzRd = MzRk / gammaM0;
  const v = Math.abs(VzEd) / VplRd;
  if (v > 0.5 && plastic) {
    const rho = (2 * v - 1) ** 2;
    const Aw = (s.hw !== undefined ? s.hw : s.h - 2 * s.tf) * s.tw;
    MyRd = Math.min(MyRd, (s.Wply - rho * Aw * Aw / (4 * s.tw)) * fy / gammaM0);
  }
  const curves = ec3BucklingCurvesI({ h: s.h, b: s.b, tf: s.tf, welded: !!s.welded, fy });
  const Ncry = Math.PI ** 2 * E * s.Iy / (Lcry * Lcry);
  const Ncrz = Math.PI ** 2 * E * s.Iz / (lz * lz);
  const lamY = Math.sqrt(NRk / Ncry);
  const lamZ = Math.sqrt(NRk / Ncrz);
  const chiY = ec3Chi(lamY, curves.y);
  const chiZ = ec3Chi(lamZ, curves.z);
  let chiLT = 1; let lamLT = 0; let Mcr = Infinity;
  if (!torsionallyRestrained) {
    Mcr = elasticCriticalMoment({ E, G, Iz: s.Iz, It: s.J, Iw: s.Iw || 0, L: llt, C1 });
    lamLT = Math.sqrt(MyRk / Mcr);
    chiLT = ec3Chi(lamLT, curves.lt);
  }
  const My = Math.abs(MyEd); const Mz = Math.abs(MzEd);
  const checks = [];
  const add = (name, uc, ref) => checks.push({ name, uc, ref });
  add('Axial resistance', Math.abs(NEd) / NplRd, 'EN 1993-1-1 §6.2.3/6.2.4');
  add('Bending about y', My / MyRd, 'EN 1993-1-1 §6.2.5 (6.2.8 with high shear)');
  add('Bending about z', Mz / MzRd, 'EN 1993-1-1 §6.2.5');
  add('Shear', v, 'EN 1993-1-1 §6.2.6');
  add('Cross-section interaction (linear sum)', Math.abs(NEd) / NplRd + My / MyRd + Mz / MzRd, 'EN 1993-1-1 §6.2.1(7)');
  add('Lateral-torsional buckling', My / (chiLT * MyRk / gammaM1), 'EN 1993-1-1 §6.3.2.2');
  const factors = { chiY, chiZ, chiLT, lambdaY: lamY, lambdaZ: lamZ, lambdaLT: lamLT, Mcr };
  if (comp > 0) {
    add('Flexural buckling about y', comp / (chiY * NRk / gammaM1), 'EN 1993-1-1 §6.3.1');
    add('Flexural buckling about z', comp / (chiZ * NRk / gammaM1), 'EN 1993-1-1 §6.3.1');
    const ny = comp / (chiY * NRk / gammaM1);
    const nz = comp / (chiZ * NRk / gammaM1);
    let kyy; let kzz; let kyz; let kzy;
    if (plastic) {
      kyy = Cmy * Math.min(1 + (lamY - 0.2) * ny, 1 + 0.8 * ny);
      kzz = Cmz * Math.min(1 + (2 * lamZ - 0.6) * nz, 1 + 1.4 * nz);
      kyz = 0.6 * kzz;
      if (torsionallyRestrained) kzy = 0.6 * kyy;
      else if (lamZ >= 0.4) kzy = Math.max(1 - 0.1 * lamZ / (CmLT - 0.25) * nz, 1 - 0.1 / (CmLT - 0.25) * nz);
      else kzy = Math.min(0.6 + lamZ, 1 - 0.1 * lamZ / (CmLT - 0.25) * nz);
    } else {
      kyy = Cmy * Math.min(1 + 0.6 * lamY * ny, 1 + 0.6 * ny);
      kzz = Cmz * Math.min(1 + 0.6 * lamZ * nz, 1 + 0.6 * nz);
      kyz = kzz;
      if (torsionallyRestrained) kzy = 0.8 * kyy;
      else kzy = Math.max(1 - 0.05 * lamZ / (CmLT - 0.25) * nz, 1 - 0.05 / (CmLT - 0.25) * nz);
    }
    const myTerm = My / (chiLT * MyRk / gammaM1);
    const mzTerm = Mz / (MzRk / gammaM1);
    add('Beam-column interaction (y)', ny + kyy * myTerm + kyz * mzTerm, 'EN 1993-1-1 Eq. (6.61), Annex B');
    add('Beam-column interaction (z)', nz + kzy * myTerm + kzz * mzTerm, 'EN 1993-1-1 Eq. (6.62), Annex B');
    Object.assign(factors, { kyy, kyz, kzy, kzz });
  }
  let gov = checks[0];
  for (const c of checks) if (c.uc > gov.uc) gov = c;
  return {
    classification: cls, curves, resistances: { NplRd, MyRd, MzRd, VplRd, NbyRd: chiY * NRk / gammaM1, NbzRd: chiZ * NRk / gammaM1, MbRd: chiLT * MyRk / gammaM1 },
    factors, checks, governing: gov, warnings,
  };
}

/* ------------------------------------------------------------------ AISC 360-16 */

/**
 * Nominal compressive strength for flexural buckling, AISC 360-16 §E3: Fe = π²E/(Lc/r)² (E3-4);
 * Fcr = 0.658^(Fy/Fe) Fy for Lc/r ≤ 4.71√(E/Fy) (E3-2), 0.877 Fe otherwise (E3-3); Pn = Fcr Ag (E3-1).
 * Non-slender elements only.
 * @param {object} o
 * @param {number} o.A gross area [m²]
 * @param {number} o.r governing radius of gyration [m]
 * @param {number} o.Lc effective length K·L [m]
 * @param {number} o.fy yield stress [Pa]
 * @param {number} [o.E=2.0e11] [Pa]
 * @returns {{Pn:number, Fcr:number, Fe:number, slenderness:number, phiPn:number, PnOmega:number}} Pn [N], stresses [Pa]; φc = 0.90, Ωc = 1.67
 */
export function aiscCompression({ A, r, Lc, fy, E = 2.0e11 }) {
  pos('A', A); pos('r', r); pos('Lc', Lc); pos('fy', fy); pos('E', E);
  const sl = Lc / r;
  const Fe = Math.PI ** 2 * E / (sl * sl);
  const Fcr = sl <= 4.71 * Math.sqrt(E / fy) ? 0.658 ** (fy / Fe) * fy : 0.877 * Fe;
  const Pn = Fcr * A;
  return { Pn, Fcr, Fe, slenderness: sl, phiPn: 0.9 * Pn, PnOmega: Pn / 1.67 };
}

/**
 * Nominal major-axis flexural strength of a doubly symmetric compact I-shape, AISC 360-16 §F2:
 * Mp = Fy Zx (F2-1); Lp = 1.76 ry √(E/Fy) (F2-5); Lr = 1.95 rts (E/(0.7Fy)) √[Jc/(Sx ho) + √((Jc/(Sx ho))² + 6.76 (0.7Fy/E)²)] (F2-6);
 * rts² = √(Iy Cw)/Sx (F2-7); c = 1; Lp < Lb ≤ Lr: Mn = Cb[Mp − (Mp − 0.7 Fy Sx)(Lb − Lp)/(Lr − Lp)] ≤ Mp (F2-2);
 * Lb > Lr: Mn = Fcr Sx ≤ Mp, Fcr = Cb π²E/(Lb/rts)² √[1 + 0.078 Jc/(Sx ho) (Lb/rts)²] (F2-3, F2-4).
 * (AISC x = strong axis = y of this library; AISC y = weak axis = z here.)
 * @param {object} o
 * @param {object} o.section iSection record
 * @param {number} o.fy [Pa]
 * @param {number} o.Lb unbraced length [m]
 * @param {number} [o.Cb=1] lateral-torsional buckling modification factor [-]
 * @param {number} [o.E=2.0e11] [Pa]
 * @returns {{Mn:number, Mp:number, Lp:number, Lr:number, rts:number, regime:string, phiMn:number, MnOmega:number}} moments [N·m], lengths [m]; φb = 0.90, Ωb = 1.67
 */
export function aiscFlexureF2({ section, fy, Lb, Cb = 1, E = 2.0e11 }) {
  const s = section;
  pos('fy', fy); pos('Cb', Cb);
  if (!(Lb >= 0)) throw new RangeError('Lb must be non-negative');
  if (!s || !(s.Wply > 0) || !(s.Iw > 0)) throw new RangeError('section must be a doubly symmetric I-section record');
  const Sx = s.Wely; const Zx = s.Wply;
  const ho = s.h - s.tf;
  const ry = s.iz;
  const rts = Math.sqrt(Math.sqrt(s.Iz * s.Iw) / Sx);
  const Mp = fy * Zx;
  const Lp = 1.76 * ry * Math.sqrt(E / fy);
  const jc = s.J / (Sx * ho);
  const Lr = 1.95 * rts * (E / (0.7 * fy)) * Math.sqrt(jc + Math.sqrt(jc * jc + 6.76 * (0.7 * fy / E) ** 2));
  let Mn; let regime;
  if (Lb <= Lp) { Mn = Mp; regime = 'yielding'; }
  else if (Lb <= Lr) { Mn = Math.min(Mp, Cb * (Mp - (Mp - 0.7 * fy * Sx) * (Lb - Lp) / (Lr - Lp))); regime = 'inelastic LTB'; }
  else {
    const q = Lb / rts;
    const Fcr = Cb * Math.PI ** 2 * E / (q * q) * Math.sqrt(1 + 0.078 * jc * q * q);
    Mn = Math.min(Mp, Fcr * Sx); regime = 'elastic LTB';
  }
  return { Mn, Mp, Lp, Lr, rts, regime, phiMn: 0.9 * Mn, MnOmega: Mn / 1.67 };
}

/**
 * AISC 360-16 §H1.1 interaction for doubly symmetric members in flexure and axial force:
 * Pr/Pc ≥ 0.2: Pr/Pc + (8/9)(Mrx/Mcx + Mry/Mcy) ≤ 1 (H1-1a); otherwise Pr/(2Pc) + (Mrx/Mcx + Mry/Mcy) ≤ 1 (H1-1b).
 * Available strengths from E3 (compression), D2 yielding (tension, φt = 0.90 / Ωt = 1.67), F2 and F6
 * (minor axis Mn = min(Fy Zy, 1.6 Fy Sy)).
 * @param {object} o
 * @param {object} o.section iSection record
 * @param {number} o.fy [Pa]
 * @param {number} [o.E=2.0e11] [Pa]
 * @param {number} [o.Pr=0] required axial strength [N], tension positive
 * @param {number} [o.Mrx=0] required major-axis flexural strength [N·m]
 * @param {number} [o.Mry=0] required minor-axis flexural strength [N·m]
 * @param {number} o.Lcx effective length for major-axis buckling [m]
 * @param {number} [o.Lcy=o.Lcx] effective length for minor-axis buckling [m]
 * @param {number} [o.Lb=o.Lcy] unbraced length for LTB [m]
 * @param {number} [o.Cb=1] [-]
 * @param {'LRFD'|'ASD'} [o.method='LRFD']
 * @returns {{Pc:number, Mcx:number, Mcy:number, uc:number, equation:string, compression:object, flexure:object}} available strengths [N, N·m], utilisation [-]
 */
export function aiscMemberCheck({ section, fy, E = 2.0e11, Pr = 0, Mrx = 0, Mry = 0, Lcx, Lcy, Lb, Cb = 1, method = 'LRFD' }) {
  const s = section;
  pos('Lcx', Lcx);
  if (method !== 'LRFD' && method !== 'ASD') throw new RangeError('method must be "LRFD" or "ASD"');
  const lcy = Lcy === undefined ? Lcx : Lcy;
  const lb = Lb === undefined ? lcy : Lb;
  const avail = (nominal, phi, omega) => (method === 'LRFD' ? phi * nominal : nominal / omega);
  const cx = aiscCompression({ A: s.A, r: s.iy, Lc: Lcx, fy, E });
  const cy = aiscCompression({ A: s.A, r: s.iz, Lc: lcy, fy, E });
  const comp = cx.Pn < cy.Pn ? cx : cy;
  const Pc = Pr < 0 ? avail(comp.Pn, 0.9, 1.67) : avail(fy * s.A, 0.9, 1.67);
  const fl = aiscFlexureF2({ section: s, fy, Lb: lb, Cb, E });
  const Mcx = avail(fl.Mn, 0.9, 1.67);
  const Mcy = avail(Math.min(fy * s.Wplz, 1.6 * fy * s.Welz), 0.9, 1.67);
  const pr = Math.abs(Pr) / Pc;
  const m = Math.abs(Mrx) / Mcx + Math.abs(Mry) / Mcy;
  const uc = pr >= 0.2 ? pr + 8 / 9 * m : pr / 2 + m;
  return { Pc, Mcx, Mcy, uc, equation: pr >= 0.2 ? 'AISC 360-16 Eq. H1-1a' : 'AISC 360-16 Eq. H1-1b', compression: comp, flexure: fl };
}

/* ------------------------------------------------------------------ serviceability */

/**
 * Deflection check against a span ratio limit (e.g. L/250 total, L/360 live — project criteria; EN 1990 A1.4 and
 * national annexes give the values).
 * @param {object} o
 * @param {number} o.deflection calculated deflection [m]
 * @param {number} o.span span (twice the length for cantilevers) [m]
 * @param {number} [o.ratio=250] limit denominator n in L/n [-]
 * @returns {{limit:number, uc:number, actualRatio:number}} limit [m], utilisation [-], L/δ [-]
 */
export function deflectionCheck({ deflection, span, ratio = 250 }) {
  pos('span', span); pos('ratio', ratio);
  if (!Number.isFinite(deflection)) throw new RangeError('deflection must be finite');
  const limit = span / ratio;
  const d = Math.abs(deflection);
  return { limit, uc: d / limit, actualRatio: d > 0 ? span / d : Infinity };
}

/* ------------------------------------------------------------------ EN 1993-1-8 connections */

/** Nominal bolt grades (EN 1993-1-8 Table 3.1): fyb, fub [Pa]. */
export const BOLT_GRADES = Object.freeze({
  4.6: { fyb: 240e6, fub: 400e6 }, 4.8: { fyb: 320e6, fub: 400e6 }, 5.6: { fyb: 300e6, fub: 500e6 }, 5.8: { fyb: 400e6, fub: 500e6 },
  6.8: { fyb: 480e6, fub: 600e6 }, 8.8: { fyb: 640e6, fub: 800e6 }, 10.9: { fyb: 900e6, fub: 1000e6 },
});

/** Tensile stress areas of coarse-thread metric bolts (ISO 898-1) [m²]. */
export const BOLT_STRESS_AREA = Object.freeze({
  M12: 84.3e-6, M16: 157e-6, M20: 245e-6, M22: 303e-6, M24: 353e-6, M27: 459e-6, M30: 561e-6, M36: 817e-6,
});

/**
 * Design resistance of a single bolt (EN 1993-1-8 Table 3.4).
 * Shear per shear plane: Fv,Rd = αv fub A/γM2 — threads in the shear plane: αv = 0.6 (grades 4.6, 5.6, 8.8) or 0.5
 * (4.8, 5.8, 6.8, 10.9) with A = As; shank in the shear plane: αv = 0.6 with the gross area.
 * Tension: Ft,Rd = k2 fub As/γM2, k2 = 0.9 (0.63 countersunk).
 * Bearing: Fb,Rd = k1 αb fu d t/γM2; αb = min(αd, fub/fu, 1) with αd = e1/(3d0) (end bolt) or p1/(3d0) − ¼ (inner);
 * k1 = min(2.8 e2/d0 − 1.7, 2.5) (edge bolt) or min(1.4 p2/d0 − 1.7, 2.5) (inner).
 * Combined shear and tension: Fv,Ed/Fv,Rd + Ft,Ed/(1.4 Ft,Rd) ≤ 1.
 * @param {object} o
 * @param {number} o.d bolt diameter [m]
 * @param {number} o.As tensile stress area [m²]
 * @param {number|string} o.grade one of BOLT_GRADES
 * @param {number} [o.shearPlanes=1]
 * @param {boolean} [o.threadsInShearPlane=true]
 * @param {number} [o.gammaM2=1.25]
 * @param {object} [o.bearing] { t plate thickness [m], fu plate ultimate strength [Pa], d0 hole diameter [m], e1, e2 [m], p1?, p2? [m] }
 * @param {number} [o.FvEd=0] design shear on the bolt [N]
 * @param {number} [o.FtEd=0] design tension on the bolt [N]
 * @returns {{FvRd:number, FtRd:number, FbRd:(number|undefined), alphaV:number, alphaB:(number|undefined), k1:(number|undefined), ucShear:number,
 *   ucTension:number, ucBearing:(number|undefined), ucCombined:number}} resistances [N], utilisations [-]
 */
export function boltResistance({ d, As, grade, shearPlanes = 1, threadsInShearPlane = true, gammaM2 = 1.25, bearing, FvEd = 0, FtEd = 0 }) {
  pos('d', d); pos('As', As); pos('gammaM2', gammaM2);
  const g = BOLT_GRADES[grade];
  if (!g) throw new RangeError(`grade "${grade}" is not a recognised bolt grade`);
  const weak = ['4.8', '5.8', '6.8', '10.9'].includes(String(grade));
  const alphaV = threadsInShearPlane ? (weak ? 0.5 : 0.6) : 0.6;
  const Ashear = threadsInShearPlane ? As : Math.PI * d * d / 4;
  const FvRd = shearPlanes * alphaV * g.fub * Ashear / gammaM2;
  const FtRd = 0.9 * g.fub * As / gammaM2;
  const out = { FvRd, FtRd, alphaV, ucShear: Math.abs(FvEd) / FvRd, ucTension: Math.abs(FtEd) / FtRd };
  out.ucCombined = Math.abs(FvEd) / FvRd + Math.abs(FtEd) / (1.4 * FtRd);
  if (bearing) {
    const { t, fu, d0, e1, e2, p1, p2 } = bearing;
    pos('bearing.t', t); pos('bearing.fu', fu); pos('bearing.d0', d0); pos('bearing.e1', e1); pos('bearing.e2', e2);
    let alphaD = e1 / (3 * d0);
    if (p1 !== undefined) alphaD = Math.min(alphaD, p1 / (3 * d0) - 0.25);
    const alphaB = Math.min(alphaD, g.fub / fu, 1);
    let k1 = Math.min(2.8 * e2 / d0 - 1.7, 2.5);
    if (p2 !== undefined) k1 = Math.min(k1, 1.4 * p2 / d0 - 1.7);
    out.alphaB = alphaB; out.k1 = k1;
    out.FbRd = k1 * alphaB * fu * d * t / gammaM2;
    out.ucBearing = Math.abs(FvEd) / out.FbRd;
  }
  return out;
}

/**
 * Elastic analysis of a bolt group under in-plane shear and torsion plus out-of-plane tension and bending
 * (classical polar-moment method: direct shear shared equally, torsional shear proportional to the distance from
 * the centroid, tension linear about the centroidal axes).
 * @param {object} o
 * @param {Array<[number, number]>} o.bolts bolt coordinates (x, y) in the faying plane [m]
 * @param {number} [o.Vx=0] shear along x through the reference point [N]
 * @param {number} [o.Vy=0] shear along y [N]
 * @param {number} [o.Mz=0] in-plane torsion about the reference point [N·m] (counter-clockwise positive)
 * @param {number} [o.N=0] tension normal to the plane [N]
 * @param {number} [o.Mx=0] bending about the x axis through the centroid [N·m] (tension at +y)
 * @param {number} [o.My=0] bending about the y axis [N·m] (tension at +x)
 * @param {[number, number]} [o.ref] reference point of the loads [m]; default the bolt-group centroid
 * @returns {{centroid:[number, number], Ip:number, forces:Array<{Fx:number, Fy:number, shear:number, tension:number}>, maxShear:number, maxTension:number}}
 */
export function boltGroupElastic({ bolts, Vx = 0, Vy = 0, Mz = 0, N = 0, Mx = 0, My = 0, ref }) {
  if (!Array.isArray(bolts) || bolts.length === 0) throw new RangeError('bolts must be a non-empty array');
  const n = bolts.length;
  let cx = 0; let cy = 0;
  for (const b of bolts) { cx += b[0]; cy += b[1]; }
  cx /= n; cy /= n;
  let Ixx = 0; let Iyy = 0;
  for (const b of bolts) { Ixx += (b[1] - cy) ** 2; Iyy += (b[0] - cx) ** 2; }
  const Ip = Ixx + Iyy;
  const r0 = ref || [cx, cy];
  // torsion about the centroid: transfer of the shear from the reference point
  const Mc = Mz + (r0[0] - cx) * Vy - (r0[1] - cy) * Vx;
  if (n > 1 && !(Ip > 0)) throw new RangeError('bolts must not all coincide');
  const forces = bolts.map((b) => {
    const dx = b[0] - cx; const dy = b[1] - cy;
    const Fx = Vx / n - (Ip > 0 ? Mc * dy / Ip : 0);
    const Fy = Vy / n + (Ip > 0 ? Mc * dx / Ip : 0);
    const tension = N / n + (Ixx > 0 ? Mx * dy / Ixx : 0) + (Iyy > 0 ? My * dx / Iyy : 0);
    return { Fx, Fy, shear: Math.hypot(Fx, Fy), tension };
  });
  return {
    centroid: [cx, cy], Ip, forces, maxShear: Math.max(...forces.map((f) => f.shear)), maxTension: Math.max(...forces.map((f) => f.tension)),
  };
}

/**
 * Correlation factor βw for fillet welds (EN 1993-1-8 Table 4.1): S235 0.80, S275 0.85, S355 0.90, S420 and S460 1.00.
 * @param {'S235'|'S275'|'S355'|'S420'|'S460'} grade
 * @returns {number} βw [-]
 */
export function weldCorrelationFactor(grade) {
  const b = { S235: 0.8, S275: 0.85, S355: 0.9, S420: 1.0, S460: 1.0 }[grade];
  if (b === undefined) throw new RangeError(`grade "${grade}" is not recognised`);
  return b;
}

/**
 * Fillet weld design resistance per unit length by the simplified method (EN 1993-1-8 §4.5.3.3):
 * fvw,d = fu/(√3 βw γM2), Fw,Rd = fvw,d · a, independent of the orientation of the force.
 * @param {object} o
 * @param {number} o.a effective throat thickness [m]
 * @param {number} o.fu nominal ultimate strength of the weaker part joined [Pa]
 * @param {number} o.betaW correlation factor [-] (see weldCorrelationFactor)
 * @param {number} [o.gammaM2=1.25]
 * @returns {{fvwd:number, FwRd:number}} design shear strength [Pa] and resistance per unit length [N/m]
 */
export function filletWeldResistance({ a, fu, betaW, gammaM2 = 1.25 }) {
  pos('a', a); pos('fu', fu); pos('betaW', betaW); pos('gammaM2', gammaM2);
  const fvwd = fu / (Math.sqrt(3) * betaW * gammaM2);
  return { fvwd, FwRd: fvwd * a };
}

/**
 * Elastic analysis of a fillet-weld group treated as lines of unit throat ("weld as a line" method): resultant
 * force per unit length at the most stressed point from direct shear, in-plane torsion, normal force and bending,
 * compared with Fw,Rd of the simplified method.
 * @param {object} o
 * @param {Array<{x1:number, y1:number, x2:number, y2:number}>} o.segments weld lines in the connection plane [m]
 * @param {number} [o.Vx=0] in-plane shear [N]
 * @param {number} [o.Vy=0] in-plane shear [N]
 * @param {number} [o.Mz=0] in-plane torsion about the weld-group centroid [N·m]
 * @param {number} [o.N=0] force normal to the plane [N]
 * @param {number} [o.Mx=0] bending about the centroidal x axis [N·m]
 * @param {number} [o.My=0] bending about the centroidal y axis [N·m]
 * @param {number} [o.FwRd] design resistance per unit length [N/m] for the utilisation
 * @returns {{length:number, centroid:[number, number], Ix:number, Iy:number, Ip:number, maxForcePerLength:number, at:[number, number], uc:(number|undefined)}}
 *   line properties (m, m³) and the governing force per unit length [N/m]
 */
export function weldGroupElastic({ segments, Vx = 0, Vy = 0, Mz = 0, N = 0, Mx = 0, My = 0, FwRd }) {
  if (!Array.isArray(segments) || segments.length === 0) throw new RangeError('segments must be a non-empty array');
  let Lt = 0; let sx = 0; let sy = 0;
  for (const s of segments) {
    const l = Math.hypot(s.x2 - s.x1, s.y2 - s.y1);
    if (!(l > 0)) throw new RangeError('weld segments must have positive length');
    Lt += l; sx += l * (s.x1 + s.x2) / 2; sy += l * (s.y1 + s.y2) / 2;
  }
  const cx = sx / Lt; const cy = sy / Lt;
  let Ix = 0; let Iy = 0;
  for (const s of segments) {
    const l = Math.hypot(s.x2 - s.x1, s.y2 - s.y1);
    const ax = s.x1 - cx; const bx = s.x2 - cx; const ay = s.y1 - cy; const by = s.y2 - cy;
    Ix += l * (ay * ay + ay * by + by * by) / 3;
    Iy += l * (ax * ax + ax * bx + bx * bx) / 3;
  }
  const Ip = Ix + Iy;
  let worst = -1; let at = [cx, cy];
  for (const s of segments) {
    for (const [x, y] of [[s.x1, s.y1], [s.x2, s.y2], [(s.x1 + s.x2) / 2, (s.y1 + s.y2) / 2]]) {
      const dx = x - cx; const dy = y - cy;
      const fx = Vx / Lt - (Ip > 0 ? Mz * dy / Ip : 0);
      const fy = Vy / Lt + (Ip > 0 ? Mz * dx / Ip : 0);
      const fz = N / Lt + (Ix > 0 ? Mx * dy / Ix : 0) + (Iy > 0 ? My * dx / Iy : 0);
      const f = Math.hypot(fx, fy, fz);
      if (f > worst) { worst = f; at = [x, y]; }
    }
  }
  return { length: Lt, centroid: [cx, cy], Ix, Iy, Ip, maxForcePerLength: worst, at, uc: FwRd === undefined ? undefined : worst / FwRd };
}

/* ------------------------------------------------------------------ lifting attachments */

/**
 * Padeye check by simple section statics (the stress formulas are those commonly used with allowable-stress lifting
 * criteria; the allowable ratios are arguments because they differ between codes — defaults are the AISC ASD /
 * API RP 2A-WSD basic values 0.9 Fy bearing, 0.4 Fy shear, 0.45 Fy tension on the net section at the pin hole,
 * 0.6 Fy combined).
 *  - bearing: σ = F/(dpin · ttot), ttot = t + 2 tc
 *  - tear-out (double shear planes from the hole to the edge): τ = F/[2 (R − dh/2) t + 2·2 (Rc − dh/2) tc]
 *  - net tension across the hole: σ = F/[(2R − dh) t + 2 (2Rc − dh) tc]
 *  - cheek-plate weld: share of the load carried by one cheek plate = F tc/ttot, over the circumferential fillet weld
 *    of throat a: τw = F tc/(ttot · π · 2Rc · a)
 *  - base section (width w, thickness t) under the in-plane load at angle from the padeye axis and an optional
 *    out-of-plane component: axial + in-plane bending + out-of-plane bending, shear, von Mises
 * F is the design sling force including all load factors (DAF, skew, consequence).
 * @param {object} o
 * @param {number} o.F design sling load [N]
 * @param {number} o.fy yield strength of the plates [Pa]
 * @param {number} o.t main plate thickness [m]
 * @param {number} o.R main plate radius about the hole centre [m]
 * @param {number} o.dHole hole diameter [m]
 * @param {number} o.dPin shackle pin diameter [m]
 * @param {number} [o.tCheek=0] thickness of each cheek plate [m]
 * @param {number} [o.rCheek=0] cheek plate radius [m]
 * @param {number} [o.aCheekWeld] throat of the cheek plate fillet weld [m]
 * @param {number} [o.baseWidth] width of the main plate at the base [m]
 * @param {number} [o.holeHeight] height of the hole centre above the base [m]
 * @param {number} [o.angleDeg=0] in-plane angle of the sling force from the padeye axis (normal to the base) [deg]
 * @param {number} [o.outOfPlaneFraction=0] lateral (out-of-plane) force as a fraction of F [-] (e.g. 0.05)
 * @param {object} [o.allow] allowable stress ratios { bearing, shear, tension, combined, weld } of fy [-]
 * @returns {{stresses:object, checks:Array<{name:string, uc:number, demand:number, capacity:number, ref:string}>, governing:object}} stresses [Pa]
 */
export function padeyeCheck({
  F, fy, t, R, dHole, dPin, tCheek = 0, rCheek = 0, aCheekWeld, baseWidth, holeHeight, angleDeg = 0, outOfPlaneFraction = 0, allow = {},
}) {
  pos('F', F); pos('fy', fy); pos('t', t); pos('R', R); pos('dHole', dHole); pos('dPin', dPin);
  if (dPin > dHole) throw new RangeError('dPin must not exceed dHole');
  if (2 * R <= dHole) throw new RangeError('R must exceed dHole/2');
  if (tCheek > 0 && !(2 * rCheek > dHole && rCheek <= R)) throw new RangeError('rCheek must lie between dHole/2 and R');
  const al = { bearing: 0.9, shear: 0.4, tension: 0.45, combined: 0.6, weld: 0.4, ...allow };
  const ttot = t + 2 * tCheek;
  const checks = [];
  const add = (name, demand, capacity, ref) => checks.push({ name, demand, capacity, uc: demand / capacity, ref });
  const bearing = F / (dPin * ttot);
  add('Pin bearing', bearing, al.bearing * fy, 'σ = F/(dpin Σt)');
  const shearArea = 2 * (R - dHole / 2) * t + (tCheek > 0 ? 4 * (rCheek - dHole / 2) * tCheek : 0);
  const tearOut = F / shearArea;
  add('Tear-out shear', tearOut, al.shear * fy, 'τ = F/(2 Σ(R − dh/2) t)');
  const netArea = (2 * R - dHole) * t + (tCheek > 0 ? 2 * (2 * rCheek - dHole) * tCheek : 0);
  const netTension = F / netArea;
  add('Net-section tension at the hole', netTension, al.tension * fy, 'σ = F/Σ(2R − dh) t');
  const stresses = { bearing, tearOut, netTension };
  if (tCheek > 0 && aCheekWeld !== undefined) {
    pos('aCheekWeld', aCheekWeld);
    const tw = F * tCheek / ttot / (2 * Math.PI * rCheek * aCheekWeld);
    stresses.cheekWeld = tw;
    add('Cheek plate weld', tw, al.weld * fy, 'τw = F tc/(Σt · 2π Rc a)');
  }
  if (baseWidth !== undefined) {
    pos('baseWidth', baseWidth); pos('holeHeight', holeHeight);
    const a = angleDeg * Math.PI / 180;
    const Fax = F * Math.cos(a); const Fsh = F * Math.sin(a); const Fout = F * outOfPlaneFraction;
    const A = baseWidth * t;
    const sAx = Fax / A;
    const sIp = Math.abs(Fsh) * holeHeight / (t * baseWidth * baseWidth / 6);
    const sOp = Fout * holeHeight / (baseWidth * t * t / 6);
    const tau = Math.hypot(Fsh, Fout) / A;
    const sig = Math.abs(sAx) + sIp + sOp;
    const vm = Math.sqrt(sig * sig + 3 * tau * tau);
    Object.assign(stresses, { baseAxial: sAx, baseInPlaneBending: sIp, baseOutOfPlaneBending: sOp, baseShear: tau, baseVonMises: vm });
    add('Base section (von Mises)', vm, al.combined * fy, 'σeq = √((σa + σb,ip + σb,op)² + 3τ²)');
  }
  let gov = checks[0];
  for (const c of checks) if (c.uc > gov.uc) gov = c;
  return { stresses, checks, governing: gov };
}

/**
 * Trunnion check: a tubular stub cantilevering from the lifted structure and loaded by the sling at a lever arm.
 * Bending σb = F e/Wel, shear τ = 2F/A (thin-wall peak) and von Mises σeq = √(σb² + 3τ²) at the root, with allowable
 * ratios supplied by the caller (defaults 0.66 Fy bending, 0.4 Fy shear, 0.75 Fy equivalent — AISC ASD / API RP 2A-WSD
 * basic values for compact tubes and combined stress).
 * @param {object} o
 * @param {number} o.F design sling load on the trunnion [N]
 * @param {number} o.arm lever arm from the sling bearing centre to the trunnion root [m]
 * @param {number} o.D trunnion outer diameter [m]
 * @param {number} o.t trunnion wall thickness [m]
 * @param {number} o.fy yield strength [Pa]
 * @param {object} [o.allow] allowable ratios { bending, shear, combined } [-]
 * @returns {{stresses:{bending:number, shear:number, vonMises:number}, checks:Array, governing:object}} stresses [Pa]
 */
export function trunnionCheck({ F, arm, D, t, fy, allow = {} }) {
  pos('F', F); pos('D', D); pos('t', t); pos('fy', fy);
  if (!(arm >= 0)) throw new RangeError('arm must be non-negative');
  if (t > D / 2) throw new RangeError('t must not exceed D/2');
  const al = { bending: 0.66, shear: 0.4, combined: 0.75, ...allow };
  const di = D - 2 * t;
  const A = Math.PI / 4 * (D * D - di * di);
  const W = Math.PI / 32 * (D ** 4 - di ** 4) / D;
  const sb = F * arm / W;
  const tau = 2 * F / A;
  const vm = Math.sqrt(sb * sb + 3 * tau * tau);
  const checks = [
    { name: 'Bending at the root', demand: sb, capacity: al.bending * fy, uc: sb / (al.bending * fy), ref: 'σb = F e/Wel' },
    { name: 'Shear', demand: tau, capacity: al.shear * fy, uc: tau / (al.shear * fy), ref: 'τ = 2F/A' },
    { name: 'Equivalent stress', demand: vm, capacity: al.combined * fy, uc: vm / (al.combined * fy), ref: 'σeq = √(σb² + 3τ²)' },
  ];
  let gov = checks[0];
  for (const c of checks) if (c.uc > gov.uc) gov = c;
  return { stresses: { bending: sb, shear: tau, vonMises: vm }, checks, governing: gov };
}
