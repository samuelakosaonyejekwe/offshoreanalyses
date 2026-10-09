/**
 * Free-span assessment in the format of DNV-RP-F105: effective span length, critical load, fundamental
 * frequency, static response, unit stress amplitudes, in-line and cross-flow VIV response models, onset
 * screening, VIV fatigue and allowable span lengths.
 *
 * Sign convention: effective axial force S_eff is tension-positive. All quantities SI.
 * The response-model construction points are those of DNV-RP-F105 (2006) §4.3–4.4, which the later
 * editions retain in form; safety factors and reduction functions are arguments so that the edition in
 * force can be applied.
 */

function pos(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}
function nonNeg(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v >= 0)) throw new RangeError(`${name} must be a non-negative finite number`);
  return v;
}
function fin(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v))) throw new RangeError(`${name} must be a finite number`);
  return v;
}

/**
 * Boundary-condition coefficients C1–C6.
 * pinned–pinned: 1.57, 1.0, 0.8, 4.93, 1/8, 5/384; fixed–fixed: 3.56, 4.0, 0.2, 14.1, 1/12, 1/384;
 * single span on seabed: 3.56, 4.0, 0.4, C4 = 14.1 (L/L_eff)² at the shoulder (8.6 at mid-span),
 * C5 = 1/(18 (L_eff/L)² − 6) at the shoulder (1/24 at mid-span), C6 = 1/384.
 *
 * Reference: DNV-RP-F105 (2006) Table 6-1.
 *
 * @param {'pinned'|'fixed'|'seabed'} type
 * @param {number} [leffOverL=1] L_eff/L for the single span on seabed [-]
 * @returns {{C1:number, C2:number, C3:number, C4:number, C4mid:number, C5:number, C5mid:number, C6:number}}
 */
export function boundaryCoefficients(type, leffOverL = 1) {
  if (type === 'pinned') return { C1: 1.57, C2: 1.0, C3: 0.8, C4: 4.93, C4mid: 4.93, C5: 1 / 8, C5mid: 1 / 8, C6: 5 / 384 };
  if (type === 'fixed') return { C1: 3.56, C2: 4.0, C3: 0.2, C4: 14.1, C4mid: 8.6, C5: 1 / 12, C5mid: 1 / 24, C6: 1 / 384 };
  if (type === 'seabed') {
    if (!(leffOverL >= 1)) throw new RangeError('leffOverL must be at least 1');
    return { C1: 3.56, C2: 4.0, C3: 0.4, C4: 14.1 / (leffOverL * leffOverL), C4mid: 8.6, C5: 1 / (18 * leffOverL * leffOverL - 6), C5mid: 1 / 24, C6: 1 / 384 };
  }
  throw new RangeError('type must be "pinned", "fixed" or "seabed"');
}

/**
 * Concrete stiffness enhancement factor CSF = k_c (EI_conc/EI_steel)^0.75.
 *
 * Reference: DNV-RP-F105 (2006) Eq. 6.5; k_c = 0.33 for asphalt and 0.25 for PP/PE corrosion coatings.
 *
 * @param {object} o
 * @param {number} o.EIconc bending stiffness of the concrete coating [N·m²]
 * @param {number} o.EIsteel bending stiffness of the steel pipe [N·m²]
 * @param {number} [o.kc=0.25] [-]
 * @returns {number} CSF [-]
 */
export function concreteStiffnessFactor({ EIconc, EIsteel, kc = 0.25 }) {
  nonNeg('EIconc', EIconc); pos('EIsteel', EIsteel); pos('kc', kc);
  return kc * (EIconc / EIsteel) ** 0.75;
}

/**
 * Effective span length ratio for a span on elastic shoulders:
 * L_eff/L = 4.73/(−0.066 β² + 1.02 β + 0.63) for β ≥ 2.7, 4.73/(0.036 β² + 0.61 β + 1.0) for β < 2.7,
 * β = log10( K L⁴/((1 + CSF) EI) ).
 *
 * Reference: DNV-RP-F105 (2006) Eq. 6.13.
 * Validity: single span on a homogeneous seabed, L/D_s below about 140.
 *
 * @param {object} o
 * @param {number} o.L span length [m]
 * @param {number} o.K soil stiffness (vertical or lateral, static or dynamic as relevant) [N/m²]
 * @param {number} o.EI steel bending stiffness [N·m²]
 * @param {number} [o.CSF=0] [-]
 * @returns {{ratio:number, Leff:number, beta:number}}
 */
export function effectiveSpanLength({ L, K, EI, CSF = 0 }) {
  pos('L', L); pos('K', K); pos('EI', EI); nonNeg('CSF', CSF);
  const beta = Math.log10((K * L ** 4) / ((1 + CSF) * EI));
  const ratio = beta >= 2.7 ? 4.73 / (-0.066 * beta * beta + 1.02 * beta + 0.63) : 4.73 / (0.036 * beta * beta + 0.61 * beta + 1.0);
  return { ratio, Leff: ratio * L, beta };
}

/**
 * Added mass coefficient near a plane boundary: C_a = 0.68 + 1.6/(1 + 5 e/D) for e/D < 0.8, else 1.
 *
 * Reference: DNV-RP-F105 (2006) Eq. 6.21.
 *
 * @param {number} gapRatio e/D, gap between pipe and seabed over outer diameter [-]
 * @returns {number} C_a [-]
 */
export function addedMassCoefficient(gapRatio) {
  nonNeg('gapRatio', gapRatio);
  return gapRatio < 0.8 ? 0.68 + 1.6 / (1 + 5 * gapRatio) : 1;
}

/**
 * Structural response of a free span: effective length, critical load, static deflection and moment,
 * fundamental frequencies in-line and cross-flow, and unit-diameter stress amplitudes.
 *
 *   P_cr = (1 + CSF) C2 π² EI/L_eff²
 *   δ = C6 q L_eff⁴/(EI (1 + CSF)) · 1/(1 + S_eff/P_cr),   M_static = C5 q L_eff²/(1 + S_eff/P_cr)
 *   f_1 = C1 √(1 + CSF) √[ EI/(m_e L_eff⁴) · (1 + S_eff/P_cr + C3 (δ/D)²) ]     (δ = 0 in-line)
 *   A_IL/CF = C4 (1 + CSF) D (D_s − t) E/L_eff²
 *
 * Reference: DNV-RP-F105 (2006) Eq. 6.8–6.13, 6.16, 6.18, Table 6-1.
 * Validity: L/D_s < 140, δ/D < 2.5, S_eff/P_cr > −0.5 for the approximate response quantities (flagged in
 * `warnings`); throws when S_eff ≤ −P_cr (the span buckles as an Euler column).
 *
 * @param {object} o
 * @param {number} o.L span length [m]
 * @param {number} o.D outer (hydrodynamic) diameter [m]
 * @param {number} o.Ds steel outer diameter [m]
 * @param {number} o.t steel wall thickness [m]
 * @param {number} o.EI steel bending stiffness [N·m²]
 * @param {number} o.me effective mass incl. contents and added mass [kg/m]
 * @param {number} [o.E=207e9] [Pa]
 * @param {number} [o.Seff=0] effective axial force, tension positive [N]
 * @param {number} [o.q=0] submerged weight (vertical static load) [N/m]
 * @param {number} [o.CSF=0] [-]
 * @param {'pinned'|'fixed'|'seabed'} [o.boundary='seabed']
 * @param {number} [o.KV] vertical soil stiffness, dynamic [N/m²] (for 'seabed', cross-flow L_eff)
 * @param {number} [o.KL] lateral soil stiffness, dynamic [N/m²] (for 'seabed', in-line L_eff; default KV)
 * @param {number} [o.KVS] static vertical soil stiffness [N/m²] (for the static L_eff; default KV)
 * @returns {{LeffCF:number, LeffIL:number, LeffStatic:number, PcrCF:number, PcrIL:number, deflection:number,
 *   momentStatic:number, stressStatic:number, f1IL:number, f1CF:number, AIL:number, ACF:number, warnings:string[]}}
 *   lengths [m]; Pcr [N]; deflection [m]; momentStatic [N·m] (shoulder value for 'seabed' and 'fixed');
 *   stressStatic bending stress at the steel outer fibre [Pa]; f1 [Hz]; A unit stress amplitudes [Pa]
 */
export function spanResponse({ L, D, Ds, t, EI, me, E = 207e9, Seff = 0, q = 0, CSF = 0, boundary = 'seabed', KV, KL, KVS }) {
  pos('L', L); pos('D', D); pos('Ds', Ds); pos('t', t); pos('EI', EI); pos('me', me); fin('Seff', Seff); nonNeg('q', q); nonNeg('CSF', CSF);
  const warnings = [];
  const ratio = (K) => (boundary === 'seabed' ? effectiveSpanLength({ L, K: pos('soil stiffness', K), EI, CSF }).ratio : 1);
  const rCF = ratio(KV);
  const rIL = ratio(KL === undefined ? KV : KL);
  const rS = ratio(KVS === undefined ? KV : KVS);
  const one = (r, withSag, delta) => {
    const c = boundaryCoefficients(boundary, r);
    const Leff = r * L;
    const Pcr = ((1 + CSF) * c.C2 * Math.PI * Math.PI * EI) / (Leff * Leff);
    const ax = 1 + Seff / Pcr;
    if (!(ax > 0)) throw new RangeError('Seff exceeds the critical load in compression: the span buckles');
    const sag = withSag ? c.C3 * (delta / D) ** 2 : 0;
    return { c, Leff, Pcr, ax, f1: c.C1 * Math.sqrt(1 + CSF) * Math.sqrt((EI / (me * Leff ** 4)) * (ax + sag)),
      A: (c.C4 * (1 + CSF) * D * (Ds - t) * E) / (Leff * Leff) };
  };
  const st = one(rS, false, 0);
  const deflection = (st.c.C6 * q * st.Leff ** 4) / (EI * (1 + CSF)) / st.ax;
  const momentStatic = (st.c.C5 * q * st.Leff ** 2) / st.ax;
  const il = one(rIL, false, 0);
  const cf = one(rCF, true, deflection);
  if (L / Ds >= 140) warnings.push('L/Ds is not below 140');
  if (deflection / D >= 2.5) warnings.push('static deflection exceeds 2.5 D');
  if (Seff / st.Pcr <= -0.5) warnings.push('Seff/Pcr is below −0.5');
  const I = EI / E;
  return { LeffCF: cf.Leff, LeffIL: il.Leff, LeffStatic: st.Leff, PcrCF: cf.Pcr, PcrIL: il.Pcr, deflection, momentStatic,
    stressStatic: (momentStatic * Ds) / (2 * I), f1IL: il.f1, f1CF: cf.f1, AIL: il.A, ACF: cf.A, warnings };
}

/**
 * Stability parameter K_s = 4π m_e ζ_T/(ρ D²) and its design value K_sd = K_s/γ_k.
 *
 * Reference: DNV-RP-F105 (2006) Eq. 4.6; total damping ζ_T = structural + soil + hydrodynamic (structural
 * 0.005 if no data, up to 0.01–0.02 with concrete coating; soil 0.005–0.02 typically; hydrodynamic 0 within
 * lock-in). γ_k = 1.0 / 1.15 / 1.30 for safety class low / normal / high (Table 2-2).
 *
 * @param {object} o
 * @param {number} o.me effective mass [kg/m]
 * @param {number} o.zetaT total modal damping ratio [-]
 * @param {number} o.D outer diameter [m]
 * @param {number} [o.rhoW=1025] [kg/m³]
 * @param {number} [o.gammaK=1.15] [-]
 * @returns {{Ks:number, Ksd:number}}
 */
export function stabilityParameter({ me, zetaT, D, rhoW = 1025, gammaK = 1.15 }) {
  pos('me', me); nonNeg('zetaT', zetaT); pos('D', D); pos('rhoW', rhoW); pos('gammaK', gammaK);
  const Ks = (4 * Math.PI * me * zetaT) / (rhoW * D * D);
  return { Ks, Ksd: Ks / gammaK };
}

function piecewise(xs, ys, x) {
  if (x <= xs[0] || x >= xs[xs.length - 1]) return 0;
  for (let i = 1; i < xs.length; i++) {
    if (x <= xs[i]) {
      const dx = xs[i] - xs[i - 1];
      return dx > 0 ? ys[i - 1] + ((ys[i] - ys[i - 1]) * (x - xs[i - 1])) / dx : Math.max(ys[i - 1], ys[i]);
    }
  }
  return 0;
}

/**
 * Turbulence-intensity / flow-angle reduction factors of the in-line response model:
 * R_Iθ,1 = 1 − π² (π/2 − √2 θ_rel)(I_c − 0.03), R_Iθ,2 = 1 − (I_c − 0.03)/0.17, both limited to [0, 1].
 *
 * Reference: DNV-RP-F105 (2006) Eq. 4.9–4.10 (as recalled; confirm before relying on the reduction — the
 * unreduced value 1.0 is conservative).
 *
 * @param {object} o
 * @param {number} o.Ic turbulence intensity [-]
 * @param {number} [o.thetaRel=Math.PI/2] angle between flow and pipe axis [rad]
 * @returns {{R1:number, R2:number}}
 */
export function inlineTurbulenceReduction({ Ic, thetaRel = Math.PI / 2 }) {
  nonNeg('Ic', Ic);
  const clamp = (v) => Math.min(1, Math.max(0, v));
  return { R1: clamp(1 - Math.PI ** 2 * (Math.PI / 2 - Math.SQRT2 * thetaRel) * (Ic - 0.03)), R2: clamp(1 - (Ic - 0.03) / 0.17) };
}

/**
 * In-line VIV response model: amplitude A_Y/D versus reduced velocity.
 * V_R,onset = 1.0/γ_on (K_sd < 0.4); (0.6 + K_sd)/γ_on (0.4 ≤ K_sd < 1.6); 2.2/γ_on (K_sd ≥ 1.6)
 * V_R,1 = 10 (A_Y,1/D) + V_R,onset;  V_R,2 = V_R,end − 2 (A_Y,2/D);  V_R,end = 4.5 − 0.8 K_sd (K_sd < 1.0), 3.7 otherwise
 * A_Y,1/D = max[ 0.18 (1 − K_sd/1.2) R_Iθ,1 ; A_Y,2/D ];  A_Y,2/D = 0.13 (1 − K_sd/1.8) R_Iθ,2
 * Curve: (V_R,onset, 0) → (V_R,1, A_Y,1) → (V_R,2, A_Y,2) → (V_R,end, 0).
 * Current-ratio factor ψ_α,IL = 0 (α < 0.5), (α − 0.5)/0.3 (0.5–0.8), 1 (α > 0.8).
 *
 * Reference: DNV-RP-F105 (2006) §4.3, Eq. 4.7–4.8, Fig. 4-1. γ_on,IL = 1.1 (Table 2-2).
 * Validity: pure in-line VIV in current-dominated flow; K_sd < 1.8.
 *
 * @param {object} o
 * @param {number} o.VR reduced velocity (U_c + U_w)/(f_n D) [-]
 * @param {number} o.Ksd design stability parameter [-]
 * @param {number} [o.gammaOn=1.1] [-]
 * @param {number} [o.R1=1] R_Iθ,1 [-]
 * @param {number} [o.R2=1] R_Iθ,2 [-]
 * @param {number} [o.alpha=1] current flow ratio U_c/(U_c + U_w) [-]
 * @returns {{AyOverD:number, psiAlpha:number, VRonset:number, VR1:number, VR2:number, VRend:number, AY1:number, AY2:number}}
 *   AyOverD is the curve value (ψ_α,IL is returned separately and enters the stress range)
 */
export function inlineResponse({ VR, Ksd, gammaOn = 1.1, R1 = 1, R2 = 1, alpha = 1 }) {
  nonNeg('VR', VR); nonNeg('Ksd', Ksd); pos('gammaOn', gammaOn);
  const VRonset = (Ksd < 0.4 ? 1.0 : Ksd < 1.6 ? 0.6 + Ksd : 2.2) / gammaOn;
  const VRend = Ksd < 1.0 ? 4.5 - 0.8 * Ksd : 3.7;
  const AY2 = Math.max(0, 0.13 * (1 - Ksd / 1.8) * R2);
  const AY1 = Math.max(0.18 * (1 - Ksd / 1.2) * R1, AY2);
  const VR1 = 10 * AY1 + VRonset;
  const VR2 = VRend - 2 * AY2;
  const psiAlpha = alpha < 0.5 ? 0 : alpha > 0.8 ? 1 : (alpha - 0.5) / 0.3;
  const AyOverD = AY1 > 0 ? piecewise([VRonset, VR1, Math.max(VR1, VR2), Math.max(VR1, VRend)], [0, AY1, AY2, 0], VR) : 0;
  return { AyOverD, psiAlpha, VRonset, VR1, VR2, VRend, AY1, AY2 };
}

/**
 * Cross-flow VIV response model: amplitude A_Z/D versus reduced velocity.
 * V_R,onset = 3 ψ_proxi ψ_mass ψ_α ψ_trench/γ_on with
 *   ψ_proxi = (4 + 1.25 e/D)/5 (e/D < 0.8), 1 otherwise;  ψ_trench = 1 + 0.5 Δ/D;
 *   ψ_mass = ½ + (ρ_s/ρ)/3 (ρ_s/ρ < 1.5), 1 otherwise;  ψ_α = 1 + α/3 (α < 0.5), 7/6 otherwise
 * V_R,1 = 7 − (7 − V_R,onset)/1.15 · (1.3 − A_Z,1/D);  V_R,2 = V_R,end − (7/1.3)(A_Z,1/D);  V_R,end = 16
 * A_Z,1/D: current-dominated (α > 0.8): 0.9 (f_ratio < 1.5), 0.9 + 0.5 (f_ratio − 1.5) (1.5–2.3), 1.3 (> 2.3);
 *          wave-dominated (α ≤ 0.8): 0.9 (KC > 30), 0.7 + 0.01 (KC − 10) (10–30), 0.7 (KC < 10);  A_Z,2 = A_Z,1
 * Curve: (2, 0) → (V_R,onset, 0.15) → (V_R,1, A_Z,1) → (V_R,2, A_Z,2) → (16, 0).
 * Damping reduction R_k = 1 − 0.15 K_sd (K_sd ≤ 4), 3.2 K_sd^−1.5 (K_sd > 4).
 *
 * Reference: DNV-RP-F105 (2006) §4.4, Eq. 4.11–4.16, Fig. 4-2. γ_on,CF = 1.2 (Table 2-2).
 * f_ratio = f_n+1,CF/f_n,CF of two consecutive cross-flow modes (single-mode response: use < 1.5).
 *
 * @param {object} o
 * @param {number} o.VR reduced velocity [-]
 * @param {number} [o.Ksd=0] [-]
 * @param {number} [o.alpha=1] current flow ratio [-]
 * @param {number} [o.KC=Infinity] Keulegan–Carpenter number U_w/(f_w D) [-]
 * @param {number} [o.fRatio=1] [-]
 * @param {number} [o.gapRatio=1] e/D [-]
 * @param {number} [o.trenchRatio=0] relative trench depth Δ/D [-]
 * @param {number} [o.massRatio=2] specific mass ρ_s/ρ [-]
 * @param {number} [o.gammaOn=1.2] [-]
 * @returns {{AzOverD:number, Rk:number, VRonset:number, VR1:number, VR2:number, VRend:number, AZ1:number}}
 *   AzOverD is the curve value before R_k
 */
export function crossflowResponse({ VR, Ksd = 0, alpha = 1, KC = Infinity, fRatio = 1, gapRatio = 1, trenchRatio = 0, massRatio = 2, gammaOn = 1.2 }) {
  nonNeg('VR', VR); nonNeg('Ksd', Ksd); nonNeg('gapRatio', gapRatio); nonNeg('trenchRatio', trenchRatio); pos('massRatio', massRatio); pos('gammaOn', gammaOn);
  const psiProxi = gapRatio < 0.8 ? (4 + 1.25 * gapRatio) / 5 : 1;
  const psiTrench = 1 + 0.5 * trenchRatio;
  const psiMass = massRatio < 1.5 ? 0.5 + massRatio / 3 : 1;
  const psiAlpha = alpha < 0.5 ? 1 + alpha / 3 : 7 / 6;
  const VRonset = (3 * psiProxi * psiMass * psiAlpha * psiTrench) / gammaOn;
  let AZ1;
  if (alpha > 0.8) AZ1 = fRatio < 1.5 ? 0.9 : fRatio <= 2.3 ? 0.9 + 0.5 * (fRatio - 1.5) : 1.3;
  else AZ1 = KC > 30 ? 0.9 : KC >= 10 ? 0.7 + 0.01 * (KC - 10) : 0.7;
  const VRend = 16;
  const VR1 = 7 - ((7 - VRonset) / 1.15) * (1.3 - AZ1);
  const VR2 = VRend - (7 / 1.3) * AZ1;
  const Rk = Ksd <= 4 ? 1 - 0.15 * Ksd : 3.2 * Ksd ** -1.5;
  const x0 = Math.min(2, VRonset);
  const AzOverD = piecewise([x0, VRonset, Math.max(VRonset, VR1), VR2, VRend], [0, 0.15, AZ1, AZ1, 0], VR);
  return { AzOverD, Rk, VRonset, VR1, VR2, VRend, AZ1 };
}

/**
 * Screening criteria for VIV fatigue.
 *   in-line:    f_n,IL/γ_IL > U_c,100yr/(V_R,onset,IL D) · (1 − (L/D)/250) · 1/ᾱ,  ᾱ = max(U_c,100/(U_w,1 + U_c,100), 0.6)
 *   cross-flow: f_n,CF/γ_CF > (U_c,100yr + U_w,1yr)/(V_R,onset,CF D)
 *
 * Reference: DNV-RP-F105 (2006) §2.3, Eq. 2.1–2.2; screening factors γ_IL = γ_CF = 1.4 (Table 2-1).
 *
 * @param {object} o
 * @param {number} o.f1IL in-line natural frequency [Hz]
 * @param {number} o.f1CF cross-flow natural frequency [Hz]
 * @param {number} o.Uc100 100-year current at pipe level [m/s]
 * @param {number} [o.Uw1=0] 1-year significant wave-induced velocity at pipe level [m/s]
 * @param {number} o.D outer diameter [m]
 * @param {number} o.L span length [m]
 * @param {number} o.VRonsetIL [-]
 * @param {number} o.VRonsetCF [-]
 * @param {number} [o.gammaIL=1.4] [-]
 * @param {number} [o.gammaCF=1.4] [-]
 * @returns {{ilOk:boolean, cfOk:boolean, ucIL:number, ucCF:number, alphaBar:number}} uc = required / available frequency
 */
export function vivScreening({ f1IL, f1CF, Uc100, Uw1 = 0, D, L, VRonsetIL, VRonsetCF, gammaIL = 1.4, gammaCF = 1.4 }) {
  pos('f1IL', f1IL); pos('f1CF', f1CF); nonNeg('Uc100', Uc100); nonNeg('Uw1', Uw1); pos('D', D); pos('L', L); pos('VRonsetIL', VRonsetIL); pos('VRonsetCF', VRonsetCF);
  const alphaBar = Uc100 + Uw1 > 0 ? Math.max(Uc100 / (Uw1 + Uc100), 0.6) : 1;
  const reqIL = (Uc100 / (VRonsetIL * D)) * Math.max(0, 1 - L / D / 250) / alphaBar;
  const reqCF = (Uc100 + Uw1) / (VRonsetCF * D);
  const ucIL = reqIL / (f1IL / gammaIL);
  const ucCF = reqCF / (f1CF / gammaCF);
  return { ilOk: ucIL < 1, cfOk: ucCF < 1, ucIL, ucCF, alphaBar };
}

function cyclesToFailure(sn, S) {
  if (typeof sn === 'function') return sn(S);
  if (sn && Number.isFinite(sn.logA) && Number.isFinite(sn.m)) return 10 ** (sn.logA - sn.m * Math.log10(S / 1e6));
  throw new RangeError('sn must be a function N(S [Pa]) or { logA, m } with S in MPa');
}

/**
 * VIV fatigue damage of a free span integrated over the long-term current distribution and a set of
 * wave conditions.
 *   in-line stress range    S_IL = 2 A_IL (A_Y/D) ψ_α,IL γ_s      at frequency f_n,IL
 *   cross-flow stress range S_CF = 2 A_CF (A_Z/D) R_k γ_s        at frequency f_n,CF
 *   damage rate = Σ_seastates p_i ∫ f_v / N(S) · dF(U_c),  life = 1/(damage per year)
 *
 * Reference: DNV-RP-F105 (2006) §2.4 and §4.3–4.4 (Eq. 4.3, 4.11); Palmgren–Miner summation. The S–N curve
 * is supplied by the caller (e.g. DNV-RP-C203), as a function N(S [Pa]) or as { logA, m } for S in MPa.
 * Cross-flow-induced in-line motion and direct wave loading (force model) are not included.
 * γ_s = 1.3, and the utilisation factor η on the damage (1.0 / 0.5 / 0.25 by safety class) is applied by
 * the caller on the returned life.
 *
 * @param {object} o
 * @param {number} o.f1IL [Hz]
 * @param {number} o.f1CF [Hz]
 * @param {number} o.AIL unit stress amplitude in-line [Pa]
 * @param {number} o.ACF unit stress amplitude cross-flow [Pa]
 * @param {number} o.D outer diameter [m]
 * @param {number} o.Ksd [-]
 * @param {{scale:number, shape:number, loc?:number}} o.current Weibull distribution of the current speed at
 *   pipe level normal to the pipe: F(U) = 1 − exp(−((U − loc)/scale)^shape) [m/s]
 * @param {Array<{prob:number, Uw:number, Tu:number}>} [o.seaStates] probability, significant wave-induced
 *   velocity at pipe level [m/s] and its period [s]; default: no waves
 * @param {function(number):number|{logA:number, m:number}} o.sn S–N curve
 * @param {number} [o.gammaS=1.3] stress range safety factor [-]
 * @param {number} [o.gammaOnIL=1.1] [-]
 * @param {number} [o.gammaOnCF=1.2] [-]
 * @param {number} [o.gapRatio=1] [-]
 * @param {number} [o.trenchRatio=0] [-]
 * @param {number} [o.massRatio=2] [-]
 * @param {number} [o.fRatio=1] [-]
 * @param {number} [o.R1=1] [-]
 * @param {number} [o.R2=1] [-]
 * @param {number} [o.nBins=200] number of current-speed bins [-]
 * @returns {{damagePerYearIL:number, damagePerYearCF:number, lifeILYears:number, lifeCFYears:number,
 *   lifeYears:number, maxStressRangeIL:number, maxStressRangeCF:number}}
 */
export function vivFatigue({ f1IL, f1CF, AIL, ACF, D, Ksd, current, seaStates, sn, gammaS = 1.3, gammaOnIL = 1.1, gammaOnCF = 1.2,
  gapRatio = 1, trenchRatio = 0, massRatio = 2, fRatio = 1, R1 = 1, R2 = 1, nBins = 200 }) {
  pos('f1IL', f1IL); pos('f1CF', f1CF); pos('AIL', AIL); pos('ACF', ACF); pos('D', D); nonNeg('Ksd', Ksd);
  if (!current) throw new RangeError('current must be given as { scale, shape, loc }');
  pos('current.scale', current.scale); pos('current.shape', current.shape);
  const loc = current.loc || 0;
  const states = seaStates && seaStates.length ? seaStates : [{ prob: 1, Uw: 0, Tu: 10 }];
  let pTot = 0;
  for (const s of states) pTot += nonNeg('seaStates.prob', s.prob);
  if (Math.abs(pTot - 1) > 1e-6) throw new RangeError('seaStates probabilities must sum to 1');
  const cdf = (u) => (u <= loc ? 0 : 1 - Math.exp(-(((u - loc) / current.scale) ** current.shape)));
  const uMax = loc + current.scale * (-Math.log(1e-9)) ** (1 / current.shape);
  const du = (uMax - loc) / nBins;
  const YEAR = 365.25 * 24 * 3600;
  let dIL = 0;
  let dCF = 0;
  let sILmax = 0;
  let sCFmax = 0;
  for (const s of states) {
    for (let i = 0; i < nBins; i++) {
      const a = loc + i * du;
      const p = cdf(a + du) - cdf(a);
      if (p <= 0) continue;
      const Uc = a + 0.5 * du;
      const U = Uc + s.Uw;
      if (!(U > 0)) continue;
      const alpha = Uc / U;
      const il = inlineResponse({ VR: U / (f1IL * D), Ksd, gammaOn: gammaOnIL, R1, R2, alpha });
      const SIL = 2 * AIL * il.AyOverD * il.psiAlpha * gammaS;
      if (SIL > 0) { dIL += (s.prob * p * f1IL * YEAR) / cyclesToFailure(sn, SIL); if (SIL > sILmax) sILmax = SIL; }
      const KC = s.Uw > 0 ? (s.Uw * s.Tu) / D : 0;
      const cf = crossflowResponse({ VR: U / (f1CF * D), Ksd, alpha, KC, fRatio, gapRatio, trenchRatio, massRatio, gammaOn: gammaOnCF });
      const SCF = 2 * ACF * cf.AzOverD * cf.Rk * gammaS;
      if (SCF > 0) { dCF += (s.prob * p * f1CF * YEAR) / cyclesToFailure(sn, SCF); if (SCF > sCFmax) sCFmax = SCF; }
    }
  }
  return { damagePerYearIL: dIL, damagePerYearCF: dCF, lifeILYears: dIL > 0 ? 1 / dIL : Infinity, lifeCFYears: dCF > 0 ? 1 / dCF : Infinity,
    lifeYears: dIL + dCF > 0 ? 1 / Math.max(dIL, dCF) : Infinity, maxStressRangeIL: sILmax, maxStressRangeCF: sCFmax };
}

/**
 * Ultimate-limit-state bending stress of a span: static (functional) bending plus the environmental bending
 * stress, the latter as half the maximum VIV stress range at the extreme velocity plus any direct wave /
 * current drag bending supplied by the caller.
 *   σ_static = M_static D_s/(2I);  σ_env = ½ max(S_IL, S_CF) + σ_drag
 *
 * Reference: DNV-RP-F105 (2006) §2.5 (peak environmental stress from the response models); the combined
 * moment should be taken into the local buckling check of the design code (strength.localBucklingLoadControlled).
 *
 * @param {object} o
 * @param {object} o.response result of spanResponse
 * @param {number} o.D [m]
 * @param {number} o.Uextreme extreme flow velocity at pipe level (current + wave) [m/s]
 * @param {number} [o.alpha=1] current flow ratio at the extreme [-]
 * @param {number} [o.Ksd=0] [-]
 * @param {number} [o.KC=Infinity] [-]
 * @param {number} [o.sigmaDrag=0] bending stress from direct drag / inertia loading [Pa]
 * @param {number} [o.gammaS=1.3] [-]
 * @param {number} [o.sigmaAllow] allowable bending stress [Pa]
 * @returns {{sigmaStatic:number, sigmaEnv:number, sigmaTotal:number, SIL:number, SCF:number, uc:(number|null)}} [Pa]
 */
export function spanUls({ response, D, Uextreme, alpha = 1, Ksd = 0, KC = Infinity, sigmaDrag = 0, gammaS = 1.3, sigmaAllow }) {
  pos('D', D); nonNeg('Uextreme', Uextreme);
  const il = inlineResponse({ VR: Uextreme / (response.f1IL * D), Ksd, alpha });
  const cf = crossflowResponse({ VR: Uextreme / (response.f1CF * D), Ksd, alpha, KC });
  const SIL = 2 * response.AIL * il.AyOverD * il.psiAlpha * gammaS;
  const SCF = 2 * response.ACF * cf.AzOverD * cf.Rk * gammaS;
  const sigmaEnv = 0.5 * Math.max(SIL, SCF) + sigmaDrag;
  const sigmaTotal = response.stressStatic + sigmaEnv;
  return { sigmaStatic: response.stressStatic, sigmaEnv, sigmaTotal, SIL, SCF, uc: sigmaAllow === undefined ? null : sigmaTotal / pos('sigmaAllow', sigmaAllow) };
}

/**
 * Maximum allowable span length (equivalently the support spacing for span rectification) for one criterion:
 *   'static'   : static bending stress σ_static(L) ≤ sigmaAllow
 *   'onset-il' : f_1,IL(L) ≥ γ_f · U/(V_R,onset,IL D)   (no in-line VIV at velocity U)
 *   'onset-cf' : f_1,CF(L) ≥ γ_f · U/(V_R,onset,CF D)   (no cross-flow VIV at velocity U)
 * solved by bisection on L (stress increases and frequency decreases monotonically with span length).
 *
 * Reference: DNV-RP-F105 (2006) §2.3 / §6.7 relations inverted numerically.
 *
 * @param {object} o
 * @param {'static'|'onset-il'|'onset-cf'} o.criterion
 * @param {object} o.span arguments of spanResponse without L
 * @param {number} [o.sigmaAllow] allowable static bending stress [Pa] ('static')
 * @param {number} [o.U] design flow velocity at pipe level [m/s] (onset criteria)
 * @param {number} [o.VRonset] onset reduced velocity [-] (onset criteria)
 * @param {number} [o.gammaF=1] safety factor on the natural frequency [-]
 * @param {number} [o.Lmin=1] [m]
 * @param {number} [o.Lmax=500] [m]
 * @returns {{L:number, response:object}} allowable span length [m] and the response at that length
 */
export function maxSpanLength({ criterion, span, sigmaAllow, U, VRonset, gammaF = 1, Lmin = 1, Lmax = 500 }) {
  let ok;
  if (criterion === 'static') {
    pos('sigmaAllow', sigmaAllow);
    ok = (r) => r.stressStatic <= sigmaAllow;
  } else if (criterion === 'onset-il' || criterion === 'onset-cf') {
    pos('U', U); pos('VRonset', VRonset);
    const fReq = (gammaF * U) / (VRonset * span.D);
    ok = criterion === 'onset-il' ? (r) => r.f1IL >= fReq : (r) => r.f1CF >= fReq;
  } else throw new RangeError('criterion must be "static", "onset-il" or "onset-cf"');
  const test = (L) => {
    try { return ok(spanResponse({ ...span, L })); } catch { return false; }
  };
  if (!test(Lmin)) throw new RangeError('the criterion is not met even at Lmin');
  let lo = Lmin;
  let hi = Lmax;
  if (test(hi)) return { L: hi, response: spanResponse({ ...span, L: hi }) };
  for (let i = 0; i < 60; i++) {
    const mid = 0.5 * (lo + hi);
    if (test(mid)) lo = mid; else hi = mid;
  }
  return { L: lo, response: spanResponse({ ...span, L: lo }) };
}
