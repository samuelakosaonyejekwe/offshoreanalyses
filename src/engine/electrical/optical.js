/**
 * Optical-fibre link engineering for submarine communication cables and fibre elements in power cables.
 *
 * References
 *  - ITU-T G.652 / G.654 / G.655 (fibre attributes), ITU-T G.650.1 (definitions), ITU-T G.957 / G.691 / G.973
 *    (link budgets, unrepeatered submarine systems), ITU-T G.977 (repeatered submarine systems).
 *  - Agrawal, Fiber-Optic Communication Systems, 4th ed. (2010): §2.2 (step-index modes, V number), §2.3–2.4
 *    (dispersion, bit-rate–distance limits), §2.5 (loss spectrum), §4.6 (Q factor and BER), §7 (OSNR).
 *  - Marcuse, "Curvature loss formula for optical fibers", J. Opt. Soc. Am. 66 (1976) 216; "Loss analysis of
 *    single-mode fiber splices", Bell Syst. Tech. J. 56 (1977) 703 (mode-field radius fit).
 *  - Smith, Appl. Opt. 11 (1972) 2489 (Raman and Brillouin thresholds).
 *  - IEC TR 62048 (power-law theory of fibre reliability).
 *
 * Unit convention: lengths in metres as everywhere in the engine, except that fibre-optic quantities which are
 * universally quoted otherwise keep their customary unit in the argument name (…DbPerKm, …Nm, …PsPerNmKm,
 * …PsPerSqrtKm, powers in dBm).
 */

import { erfc, besselJ0, besselJ1, besselK0, besselK1 } from '../math/special.js';

const C_LIGHT = 299792458;

function pos(v, name) {
  if (!(v > 0) || !Number.isFinite(v)) throw new RangeError(`${name}: must be a positive finite number`);
}
function nonNeg(v, name) {
  if (!(v >= 0) || !Number.isFinite(v)) throw new RangeError(`${name}: must be a non-negative finite number`);
}

/**
 * Typical cabled-fibre attenuation by ITU-T fibre category and wavelength, dB/km. These are representative
 * modern values (not limits of the Recommendations) and are meant to be edited: G.652.D 0.35 / 0.20 / 0.23 at
 * 1310 / 1550 / 1625 nm; G.654 (cut-off shifted, pure-silica core, for long unrepeatered and submarine spans)
 * 0.17 at 1550 nm; G.655 (non-zero dispersion shifted) 0.22 at 1550 nm. Missing entries are null.
 */
export const FIBRE_ATTENUATION_DB_PER_KM = Object.freeze({
  'G.652': Object.freeze({ 1310: 0.35, 1550: 0.20, 1625: 0.23 }),
  'G.654': Object.freeze({ 1310: null, 1550: 0.17, 1625: null }),
  'G.655': Object.freeze({ 1310: null, 1550: 0.22, 1625: null }),
});

/**
 * Optical power budget of a link.
 * Total loss = α·L + nSplices·splice + nConnectors·connector + repair margin + ageing + other; received power =
 * launch power + Raman/ROPA gain − total loss; margin = received − sensitivity − system margin.
 * The maximum span is the length at which the margin is zero, with splice count growing with length when a
 * splice interval (cable section length) is given.
 * Reference: ITU-T G.957 Appendix / G.973 budget method.
 * @param {object} o
 * @param {number} o.length link length, m
 * @param {number} o.attenuationDbPerKm fibre attenuation at the operating wavelength, dB/km
 * @param {number} o.launchPowerDbm transmitter output, dBm
 * @param {number} o.sensitivityDbm receiver sensitivity at the required BER, dBm
 * @param {number} [o.spliceLossDb=0.05] loss per splice, dB @param {number} [o.nSplices] number of splices
 * @param {number} [o.spliceInterval] distance between splices (factory/joint length), m (used when nSplices is
 *        omitted: nSplices = ceil(L/interval) − 1 + 2 terminal splices)
 * @param {number} [o.connectorLossDb=0.5] @param {number} [o.nConnectors=2]
 * @param {number} [o.repairMarginDb=0] allowance for future repairs (extra cable and splices), dB
 * @param {number} [o.ageingDb=0] fibre and component ageing (e.g. hydrogen), dB
 * @param {number} [o.systemMarginDb=3] operating margin kept in hand, dB
 * @param {number} [o.otherLossDb=0] e.g. dispersion penalty, WDM filters, dB
 * @param {number} [o.amplifierGainDb=0] distributed Raman / remote optically pumped amplifier budget improvement, dB
 * @returns {{fibreLossDb:number, spliceLossDb:number, connectorLossDb:number, totalLossDb:number,
 *            receivedPowerDbm:number, marginDb:number, ok:boolean, nSplices:number, maxLength:number}}
 *          maxLength in m
 */
export function linkBudget({ length, attenuationDbPerKm, launchPowerDbm, sensitivityDbm, spliceLossDb = 0.05, nSplices,
  spliceInterval, connectorLossDb = 0.5, nConnectors = 2, repairMarginDb = 0, ageingDb = 0, systemMarginDb = 3,
  otherLossDb = 0, amplifierGainDb = 0 }) {
  nonNeg(length, 'length'); pos(attenuationDbPerKm, 'attenuationDbPerKm');
  if (!Number.isFinite(launchPowerDbm) || !Number.isFinite(sensitivityDbm)) throw new RangeError('launchPowerDbm, sensitivityDbm: must be finite');
  nonNeg(spliceLossDb, 'spliceLossDb'); nonNeg(connectorLossDb, 'connectorLossDb'); nonNeg(nConnectors, 'nConnectors');
  if (spliceInterval !== undefined) pos(spliceInterval, 'spliceInterval');
  const splicesFor = (L) => {
    if (nSplices !== undefined) return nSplices;
    if (spliceInterval === undefined) return 0;
    return Math.max(0, Math.ceil(L / spliceInterval - 1e-9) - 1) + 2;
  };
  const fixed = nConnectors * connectorLossDb + repairMarginDb + ageingDb + otherLossDb;
  const lossAt = (L) => (attenuationDbPerKm * L) / 1000 + splicesFor(L) * spliceLossDb + fixed;
  const marginAt = (L) => launchPowerDbm + amplifierGainDb - lossAt(L) - sensitivityDbm - systemMarginDb;
  const fibreLossDb = (attenuationDbPerKm * length) / 1000;
  const ns = splicesFor(length);
  const totalLossDb = lossAt(length);
  // Maximum span: margin is non-increasing in L; bisection copes with the splice-count steps.
  let maxLength = 0;
  if (marginAt(0) > 0) {
    let lo = 0; let hi = 1000;
    while (marginAt(hi) > 0 && hi < 1e9) hi *= 2;
    for (let i = 0; i < 100 && hi - lo > 1e-3; i++) { const mid = 0.5 * (lo + hi); if (marginAt(mid) >= 0) lo = mid; else hi = mid; }
    maxLength = lo;
  }
  const receivedPowerDbm = launchPowerDbm + amplifierGainDb - totalLossDb;
  const marginDb = receivedPowerDbm - sensitivityDbm - systemMarginDb;
  return { fibreLossDb, spliceLossDb: ns * spliceLossDb, connectorLossDb: nConnectors * connectorLossDb, totalLossDb,
    receivedPowerDbm, marginDb, ok: marginDb >= 0, nSplices: ns, maxLength };
}

/**
 * Spectral attenuation of silica fibre: α(λ) = A/λ⁴ + B + C·exp(−D/λ) + OH peak, λ in µm, α in dB/km.
 * Rayleigh scattering A/λ⁴ (A ≈ 0.75–1.0 dB·µm⁴/km for germanium-doped cores, lower for pure-silica cores),
 * wavelength-independent imperfection loss B, infrared absorption tail C·exp(−D/λ) with the commonly used fit
 * C = 7.81e11 dB/km, D = 48.48 µm, and a Gaussian OH overtone peak at 1.383 µm (zero for low-water-peak fibre).
 * All coefficients are arguments; the defaults are representative, not a specification.
 * Reference: Agrawal §2.5; ITU-T G.652 Appendix (spectral model form).
 * @param {object} o
 * @param {number} o.wavelengthNm wavelength, nm (about 1000–1700)
 * @param {number} [o.rayleigh=0.9] A, dB·µm⁴/km @param {number} [o.constant=0.01] B, dB/km
 * @param {number} [o.irAmplitude=7.81e11] C, dB/km @param {number} [o.irExponent=48.48] D, µm
 * @param {number} [o.ohPeak=0] OH peak height at 1383 nm, dB/km @param {number} [o.ohWidthNm=15] Gaussian σ, nm
 * @returns {{total:number, rayleigh:number, infrared:number, oh:number}} dB/km
 */
export function spectralAttenuation({ wavelengthNm, rayleigh = 0.9, constant = 0.01, irAmplitude = 7.81e11, irExponent = 48.48,
  ohPeak = 0, ohWidthNm = 15 }) {
  pos(wavelengthNm, 'wavelengthNm'); nonNeg(rayleigh, 'rayleigh'); nonNeg(constant, 'constant'); nonNeg(ohPeak, 'ohPeak'); pos(ohWidthNm, 'ohWidthNm');
  const um = wavelengthNm / 1000;
  const ray = rayleigh / um ** 4;
  const ir = irAmplitude * Math.exp(-irExponent / um);
  const oh = ohPeak * Math.exp(-0.5 * ((wavelengthNm - 1383) / ohWidthNm) ** 2);
  return { total: ray + constant + ir + oh, rayleigh: ray, infrared: ir, oh };
}

/**
 * Chromatic dispersion of standard single-mode fibre from the three-term Sellmeier fit:
 * D(λ) = S0/4·(λ − λ0⁴/λ³) ps/(nm·km), and the group-velocity dispersion β2 = −D·λ²/(2πc).
 * Reference: ITU-T G.652 / G.650.1 (1260–1360 nm fit, commonly extrapolated to 1550 nm, where it gives about
 * 17 ps/(nm·km) for λ0 = 1310 nm and S0 = 0.092 ps/(nm²·km)).
 * @param {object} o
 * @param {number} o.wavelengthNm nm
 * @param {number} [o.zeroDispersionNm=1310] λ0 (G.652: 1300–1324 nm)
 * @param {number} [o.slopePsPerNm2Km=0.092] S0 at λ0 (G.652 maximum 0.092; typical 0.086)
 * @returns {{DPsPerNmKm:number, beta2:number}} dispersion in ps/(nm·km) and β2 in s²/m
 */
export function chromaticDispersion({ wavelengthNm, zeroDispersionNm = 1310, slopePsPerNm2Km = 0.092 }) {
  pos(wavelengthNm, 'wavelengthNm'); pos(zeroDispersionNm, 'zeroDispersionNm'); nonNeg(slopePsPerNm2Km, 'slopePsPerNm2Km');
  const D = (slopePsPerNm2Km / 4) * (wavelengthNm - zeroDispersionNm ** 4 / wavelengthNm ** 3);
  // ps/(nm·km) = 1e-6 s/m²
  const lam = wavelengthNm * 1e-9;
  return { DPsPerNmKm: D, beta2: (-(D * 1e-6) * lam * lam) / (2 * Math.PI * C_LIGHT) };
}

/**
 * Dispersion-limited length.
 *  - Source with spectral width much larger than the modulation bandwidth (rms width σλ): B·L·|D|·σλ ≤ 1/4.
 *  - Narrow-linewidth source, chirp-free externally modulated NRZ (transform-limited pulses): B²·|β2|·L ≤ 1/16.
 * Reference: Agrawal §2.4.3 (Eqs. for the large and small spectral-width limits).
 * @param {object} o
 * @param {number} o.bitRate bit rate B, bit/s
 * @param {number} o.DPsPerNmKm dispersion, ps/(nm·km)
 * @param {number} o.wavelengthNm nm
 * @param {number} [o.sourceWidthNm] rms spectral width σλ of a broad source, nm (omit for a narrow source)
 * @returns {{length:number, criterion:string, beta2:number}} m (Infinity for zero dispersion), s²/m
 */
export function dispersionLimitedLength({ bitRate, DPsPerNmKm, wavelengthNm, sourceWidthNm }) {
  pos(bitRate, 'bitRate'); pos(wavelengthNm, 'wavelengthNm');
  if (!Number.isFinite(DPsPerNmKm)) throw new RangeError('DPsPerNmKm: must be finite');
  const lam = wavelengthNm * 1e-9;
  const Dsi = Math.abs(DPsPerNmKm) * 1e-6; // s/m²
  const beta2 = (-(DPsPerNmKm * 1e-6) * lam * lam) / (2 * Math.PI * C_LIGHT);
  if (Dsi === 0) return { length: Infinity, criterion: 'none', beta2 };
  if (sourceWidthNm !== undefined) {
    pos(sourceWidthNm, 'sourceWidthNm');
    return { length: 1 / (4 * bitRate * Dsi * sourceWidthNm * 1e-9), criterion: 'B·L·|D|·σλ ≤ 1/4', beta2 };
  }
  return { length: 1 / (16 * bitRate * bitRate * Math.abs(beta2)), criterion: 'B²·|β2|·L ≤ 1/16', beta2 };
}

/**
 * PMD-limited length: mean differential group delay DGD = PMDcoef·sqrt(L) must stay below a fraction of the bit
 * period (0.1 by default, the usual NRZ allowance for about 1 dB penalty).
 * Reference: ITU-T G.691 Appendix I; Agrawal §2.4.5.
 * @param {object} o
 * @param {number} o.bitRate bit/s @param {number} o.pmdPsPerSqrtKm PMD coefficient, ps/√km (G.652.D ≤ 0.20)
 * @param {number} [o.fraction=0.1] allowed DGD as a fraction of the bit period
 * @param {number} [o.length] m, to report the DGD of an actual link
 * @returns {{length:number, dgdPs:?number}} limiting length (m) and DGD of the given link (ps)
 */
export function pmdLimitedLength({ bitRate, pmdPsPerSqrtKm, fraction = 0.1, length }) {
  pos(bitRate, 'bitRate'); pos(pmdPsPerSqrtKm, 'pmdPsPerSqrtKm'); pos(fraction, 'fraction');
  const Tps = 1e12 / bitRate;
  const km = ((fraction * Tps) / pmdPsPerSqrtKm) ** 2;
  return { length: km * 1000, dgdPs: length !== undefined ? pmdPsPerSqrtKm * Math.sqrt(length / 1000) : null };
}

/**
 * Length of dispersion-compensating fibre that cancels the accumulated dispersion of a span:
 * L_dcf = −D·L/D_dcf, and the residual dispersion for a given compensation ratio.
 * @param {object} o
 * @param {number} o.length span length, m @param {number} o.DPsPerNmKm span fibre dispersion, ps/(nm·km)
 * @param {number} o.DdcfPsPerNmKm dispersion of the compensating fibre (negative, typically −80 … −150)
 * @param {number} [o.ratio=1] fraction of the accumulated dispersion to compensate
 * @returns {{lengthDcf:number, accumulatedPsPerNm:number, residualPsPerNm:number}} m, ps/nm, ps/nm
 */
export function dispersionCompensation({ length, DPsPerNmKm, DdcfPsPerNmKm, ratio = 1 }) {
  nonNeg(length, 'length'); nonNeg(ratio, 'ratio');
  if (!(DdcfPsPerNmKm * DPsPerNmKm < 0)) throw new RangeError('DdcfPsPerNmKm: must have the opposite sign to DPsPerNmKm');
  const acc = (DPsPerNmKm * length) / 1000;
  return { lengthDcf: ((-ratio * acc) / DdcfPsPerNmKm) * 1000, accumulatedPsPerNm: acc, residualPsPerNm: acc * (1 - ratio) };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Waveguide theory                                                                                             */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Step-index fibre parameters: NA = sqrt(n1² − n2²), V = 2π·a·NA/λ, single-mode cut-off wavelength
 * λc = 2π·a·NA/2.405 (first zero of J0), relative index difference Δ = (n1² − n2²)/(2·n1²), approximate number
 * of guided modes V²/2 for V ≫ 1 (two polarisations counted; 2 for single-mode operation), and the Gaussian
 * mode-field diameter from Marcuse's fit w/a = 0.65 + 1.619·V^(−3/2) + 2.879·V^(−6) (0.8 < V < 2.5).
 * Reference: Agrawal §2.2; Marcuse (1977).
 * @param {object} o
 * @param {number} o.coreRadius a, m @param {number} o.n1 core index @param {number} o.n2 cladding index
 * @param {number} o.wavelengthNm nm
 * @returns {{NA:number, V:number, cutoffWavelengthNm:number, delta:number, singleMode:boolean, modes:number,
 *            modeFieldDiameter:number, acceptanceAngle:number}} MFD in m (NaN outside the fit range), acceptance
 *          half-angle in air in rad
 */
export function stepIndexFibre({ coreRadius, n1, n2, wavelengthNm }) {
  pos(coreRadius, 'coreRadius'); pos(n1, 'n1'); pos(n2, 'n2'); pos(wavelengthNm, 'wavelengthNm');
  if (!(n1 > n2)) throw new RangeError('n1: core index must exceed cladding index n2');
  const NA = Math.sqrt(n1 * n1 - n2 * n2);
  const lam = wavelengthNm * 1e-9;
  const V = (2 * Math.PI * coreRadius * NA) / lam;
  const singleMode = V < 2.405;
  const w = coreRadius * (0.65 + 1.619 * Math.pow(V, -1.5) + 2.879 * Math.pow(V, -6));
  return { NA, V, cutoffWavelengthNm: ((2 * Math.PI * coreRadius * NA) / 2.405) * 1e9, delta: (n1 * n1 - n2 * n2) / (2 * n1 * n1),
    singleMode, modes: singleMode ? 2 : Math.round((V * V) / 2), modeFieldDiameter: V > 0.8 && V < 2.5 ? 2 * w : NaN,
    acceptanceAngle: Math.asin(Math.min(1, NA)) };
}

/**
 * Power-law (α-profile) refractive index: n(r) = n1·sqrt(1 − 2Δ·(r/a)^g) for r ≤ a, n2 = n1·sqrt(1 − 2Δ) beyond;
 * g → ∞ is the step profile, g = 2 the parabolic profile. Guided-mode count ≈ g/(g + 2)·V²/2.
 * Reference: Gloge & Marcatili, Bell Syst. Tech. J. 52 (1973) 1563.
 * @param {object} o
 * @param {number} o.r radius, m @param {number} o.coreRadius a, m @param {number} o.n1 index on the axis
 * @param {number} o.delta relative index difference Δ @param {number} [o.g=Infinity] profile exponent
 * @returns {number} refractive index at r
 */
export function refractiveIndexProfile({ r, coreRadius, n1, delta, g = Infinity }) {
  nonNeg(r, 'r'); pos(coreRadius, 'coreRadius'); pos(n1, 'n1');
  // Infinity is the documented default (step profile), so it must pass the positivity check
  if (g !== Infinity) pos(g, 'g');
  if (!(delta > 0 && delta < 0.5)) throw new RangeError('delta: must be in (0, 0.5)');
  if (r >= coreRadius) return n1 * Math.sqrt(1 - 2 * delta);
  const x = g === Infinity ? 0 : Math.pow(r / coreRadius, g);
  return n1 * Math.sqrt(1 - 2 * delta * x);
}

/**
 * LP01 mode of a weakly guiding step-index fibre: solves u·J1(u)/J0(u) = w·K1(w)/K0(w) with u² + w² = V².
 * Reference: Gloge, Appl. Opt. 10 (1971) 2252.
 * @param {number} V normalised frequency (> 0)
 * @returns {{u:number, w:number, b:number}} transverse parameters and normalised propagation constant b = w²/V²
 */
export function lp01Mode(V) {
  pos(V, 'V');
  const f = (u) => {
    const w = Math.sqrt(V * V - u * u);
    return u * besselJ1(u) * besselK0(w) - w * besselK1(w) * besselJ0(u);
  };
  // The LP01 root lies below both V and the first zero of J0.
  let lo = 1e-9; let hi = Math.min(V, 2.404825557695773) * (1 - 1e-12);
  let flo = f(lo);
  for (let i = 0; i < 200; i++) {
    const mid = 0.5 * (lo + hi);
    const fm = f(mid);
    if ((fm < 0) === (flo < 0)) { lo = mid; flo = fm; } else hi = mid;
  }
  const u = 0.5 * (lo + hi);
  const w = Math.sqrt(V * V - u * u);
  return { u, w, b: (w * w) / (V * V) };
}

/**
 * Macro-bend loss of the LP01 mode of a step-index fibre (Marcuse 1976):
 *   2α = sqrt(π)·κ²·exp(−(2/3)·(γ³/β²)·R) / (2·γ^(3/2)·V²·sqrt(R)·K1(γa)²)   [1/m, power],
 * with κ = u/a, γ = w/a, β the propagation constant, R the bend radius; loss in dB/m = 4.343·2α.
 * The formula treats an infinite cladding; real coated fibres show oscillations about this value and
 * bend-insensitive (trench-assisted) designs are far better, so use it for trends and for the minimum bend
 * radius of the cable, and manufacturer macro-bend data (ITU-T G.652/G.657 test radii) for budgets.
 * Reference: Marcuse (1976), J. Opt. Soc. Am. 66, 216 (general LPνm formula evaluated for ν = 0, eν = 2).
 * @param {object} o
 * @param {number} o.bendRadius R, m @param {number} o.coreRadius a, m @param {number} o.n1 @param {number} o.n2
 * @param {number} o.wavelengthNm nm @param {number} [o.turns] number of full turns for the total loss
 * @returns {{lossDbPerM:number, lossDbPerTurn:number, lossDb:?number, V:number}}
 */
export function macroBendLoss({ bendRadius, coreRadius, n1, n2, wavelengthNm, turns }) {
  pos(bendRadius, 'bendRadius'); pos(coreRadius, 'coreRadius'); pos(wavelengthNm, 'wavelengthNm');
  if (!(n1 > n2 && n2 > 0)) throw new RangeError('n1: core index must exceed cladding index n2');
  const lam = wavelengthNm * 1e-9;
  const k0 = (2 * Math.PI) / lam;
  const V = k0 * coreRadius * Math.sqrt(n1 * n1 - n2 * n2);
  const { u, w } = lp01Mode(V);
  const kappa = u / coreRadius; const gamma = w / coreRadius;
  const beta = Math.sqrt(n1 * n1 * k0 * k0 - kappa * kappa);
  const k1 = besselK1(w);
  const twoAlpha = (Math.sqrt(Math.PI) * kappa * kappa * Math.exp((-2 / 3) * (gamma ** 3 / (beta * beta)) * bendRadius))
    / (2 * Math.pow(gamma, 1.5) * V * V * Math.sqrt(bendRadius) * k1 * k1);
  const lossDbPerM = 10 * Math.LOG10E * twoAlpha;
  const lossDbPerTurn = lossDbPerM * 2 * Math.PI * bendRadius;
  return { lossDbPerM, lossDbPerTurn, lossDb: turns !== undefined ? lossDbPerTurn * turns : null, V };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Amplified systems                                                                                            */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * OSNR of a chain of identical amplified spans (0.1 nm reference bandwidth at 1550 nm):
 * OSNR[dB] = 58 + Pch[dBm] − NF[dB] − Lspan[dB] − 10·log10(Nspans).
 * The constant 58 dB is −10·log10(h·ν·Δν/1 mW) for Δν = 12.5 GHz at 193.4 THz.
 * Reference: ITU-T G.692 Appendix / G.977; Agrawal §7.
 * @param {object} o
 * @param {number} o.channelPowerDbm launch power per channel, dBm @param {number} o.noiseFigureDb amplifier NF, dB
 * @param {number} o.spanLossDb loss of one span, dB @param {number} o.nSpans number of spans
 * @returns {number} OSNR, dB
 */
export function osnr({ channelPowerDbm, noiseFigureDb, spanLossDb, nSpans }) {
  nonNeg(noiseFigureDb, 'noiseFigureDb'); nonNeg(spanLossDb, 'spanLossDb');
  if (!(nSpans >= 1)) throw new RangeError('nSpans: must be ≥ 1');
  return 58 + channelPowerDbm - noiseFigureDb - spanLossDb - 10 * Math.log10(nSpans);
}

/**
 * Smallest number of equal amplified spans (repeaters) that meets a required OSNR over a route.
 * Span loss = α·L/N + extra loss per span. Returns the repeater count N − 1 (terminal amplifiers excluded) and
 * the repeater spacing.
 * @param {object} o
 * @param {number} o.length route length, m @param {number} o.attenuationDbPerKm dB/km
 * @param {number} o.channelPowerDbm dBm @param {number} o.noiseFigureDb dB
 * @param {number} o.requiredOsnrDb required OSNR at the receiver including margins, dB
 * @param {number} [o.extraSpanLossDb=0] splices, repair and ageing allowance per span, dB
 * @param {number} [o.maxSpans=2000]
 * @returns {{nSpans:number, nRepeaters:number, spacing:number, spanLossDb:number, osnrDb:number}} spacing in m
 */
export function repeaterSpacing({ length, attenuationDbPerKm, channelPowerDbm, noiseFigureDb, requiredOsnrDb, extraSpanLossDb = 0, maxSpans = 2000 }) {
  pos(length, 'length'); pos(attenuationDbPerKm, 'attenuationDbPerKm'); nonNeg(extraSpanLossDb, 'extraSpanLossDb');
  for (let n = 1; n <= maxSpans; n++) {
    const spanLossDb = (attenuationDbPerKm * length) / 1000 / n + extraSpanLossDb;
    const o = osnr({ channelPowerDbm, noiseFigureDb, spanLossDb, nSpans: n });
    if (o >= requiredOsnrDb) return { nSpans: n, nRepeaters: n - 1, spacing: length / n, spanLossDb, osnrDb: o };
  }
  throw new RangeError('requiredOsnrDb: cannot be met within maxSpans (raise launch power or lower the requirement)');
}

/**
 * Q factor ↔ bit-error ratio for binary detection with Gaussian noise: BER = ½·erfc(Q/√2); Q[dB] = 20·log10(Q).
 * Give either Q (linear), QDb or BER; the others are returned. A FEC threshold (pre-FEC BER limit, e.g. 3.8e-3
 * for 7 % hard-decision FEC, about 2e-2 for soft-decision FEC — code-specific, argument) gives the Q margin.
 * Reference: Agrawal §4.6.
 * @param {object} o
 * @param {number} [o.Q] @param {number} [o.QDb] @param {number} [o.ber]
 * @param {number} [o.fecLimitBer] pre-FEC BER threshold of the code in use
 * @returns {{Q:number, QDb:number, ber:number, fecLimitQDb:?number, marginDb:?number}}
 */
export function qFactorBer({ Q, QDb, ber, fecLimitBer }) {
  const berOf = (q) => 0.5 * erfc(q / Math.SQRT2);
  const qOf = (b) => {
    if (!(b > 0 && b < 0.5)) throw new RangeError('ber: must be in (0, 0.5)');
    let lo = 0; let hi = 40;
    for (let i = 0; i < 200; i++) { const mid = 0.5 * (lo + hi); if (berOf(mid) > b) lo = mid; else hi = mid; }
    return 0.5 * (lo + hi);
  };
  let q;
  if (Q !== undefined) { pos(Q, 'Q'); q = Q; } else if (QDb !== undefined) q = Math.pow(10, QDb / 20);
  else if (ber !== undefined) q = qOf(ber);
  else throw new RangeError('Q, QDb or ber: one of them is required');
  const out = { Q: q, QDb: 20 * Math.log10(q), ber: ber !== undefined && Q === undefined && QDb === undefined ? ber : berOf(q),
    fecLimitQDb: null, marginDb: null };
  if (fecLimitBer !== undefined) {
    out.fecLimitQDb = 20 * Math.log10(qOf(fecLimitBer));
    out.marginDb = out.QDb - out.fecLimitQDb;
  }
  return out;
}

/**
 * Effective interaction length and threshold powers of stimulated Brillouin and Raman scattering:
 * Leff = (1 − e^(−αL))/α; P_SBS ≈ 21·Aeff/(gB·Leff), P_SRS ≈ 16·Aeff/(gR·Leff) (single unmodulated channel).
 * Gains are arguments (silica near 1550 nm: gB ≈ 5e-11 m/W for a narrow-linewidth source, gR ≈ 6e-14–1e-13 m/W);
 * the launch power per channel is normally kept well below these and is in practice limited by Kerr effects
 * (self- and cross-phase modulation), which need a propagation simulation.
 * Reference: Smith (1972); Agrawal §2.6.
 * @param {object} o
 * @param {number} o.length m @param {number} o.attenuationDbPerKm dB/km
 * @param {number} o.effectiveArea Aeff, m² (about 80e-12 for G.652, 110–150e-12 for G.654)
 * @param {number} [o.brillouinGain=5e-11] m/W @param {number} [o.ramanGain=1e-13] m/W
 * @returns {{effectiveLength:number, sbsThresholdW:number, srsThresholdW:number, sbsThresholdDbm:number,
 *            srsThresholdDbm:number}}
 */
export function nonlinearThresholds({ length, attenuationDbPerKm, effectiveArea, brillouinGain = 5e-11, ramanGain = 1e-13 }) {
  pos(length, 'length'); pos(attenuationDbPerKm, 'attenuationDbPerKm'); pos(effectiveArea, 'effectiveArea');
  pos(brillouinGain, 'brillouinGain'); pos(ramanGain, 'ramanGain');
  const alpha = (attenuationDbPerKm / 1000) / (10 * Math.LOG10E); // 1/m
  const effectiveLength = -Math.expm1(-alpha * length) / alpha;
  const sbs = (21 * effectiveArea) / (brillouinGain * effectiveLength);
  const srs = (16 * effectiveArea) / (ramanGain * effectiveLength);
  const dbm = (p) => 10 * Math.log10(p / 1e-3);
  return { effectiveLength, sbsThresholdW: sbs, srsThresholdW: srs, sbsThresholdDbm: dbm(sbs), srsThresholdDbm: dbm(srs) };
}

/**
 * Power-feed equipment (PFE) voltage of a repeatered system fed with constant current through the cable's
 * power conductor: V = N·Vrepeater + I·r·L + branching-unit drops + earth-potential allowance.
 * Double-end feeding shares the voltage between the two stations (+V/2 and −V/2), each normally rated so that
 * one station can feed the whole system alone for a limited period (that requirement is a project decision).
 * Reference: ITU-T G.977 §power feeding.
 * @param {object} o
 * @param {number} o.nRepeaters @param {number} o.repeaterVoltage voltage drop per repeater, V
 * @param {number} o.current line current, A (typically 0.6–1.6 A)
 * @param {number} o.conductorResistance DC resistance of the power conductor, Ω/m (about 0.7–1.6 Ω/km)
 * @param {number} o.length system length, m
 * @param {number} [o.earthPotentialMargin=0] allowance for earth potential differences (order 0.1 V/km), V
 * @param {number} [o.otherDrops=0] branching units, equalisers, V
 * @param {'single'|'double'} [o.feed='double']
 * @returns {{totalVoltage:number, voltagePerStation:number, cableDrop:number, repeaterDrop:number, power:number}}
 *          V, V, V, V, W
 */
export function powerFeedVoltage({ nRepeaters, repeaterVoltage, current, conductorResistance, length, earthPotentialMargin = 0,
  otherDrops = 0, feed = 'double' }) {
  nonNeg(nRepeaters, 'nRepeaters'); nonNeg(repeaterVoltage, 'repeaterVoltage'); pos(current, 'current');
  nonNeg(conductorResistance, 'conductorResistance'); nonNeg(length, 'length');
  if (feed !== 'single' && feed !== 'double') throw new RangeError('feed: must be single or double');
  const repeaterDrop = nRepeaters * repeaterVoltage;
  const cableDrop = current * conductorResistance * length;
  const totalVoltage = repeaterDrop + cableDrop + otherDrops + earthPotentialMargin;
  return { totalVoltage, voltagePerStation: feed === 'double' ? totalVoltage / 2 : totalVoltage, cableDrop, repeaterDrop,
    power: totalVoltage * current };
}

/**
 * Reach of a distributed fibre sensing interrogator (DTS, DAS, BOTDR/BOTDA) from its optical dynamic range:
 * range = dynamic range / (α·passes), with passes = 2 for single-ended (backscatter) schemes and 1 for loop
 * (double-ended) schemes measured per direction. Splice and connector losses reduce the available range.
 * A hook for the monitoring analyses: the interrogator dynamic range is equipment data.
 * @param {object} o
 * @param {number} o.dynamicRangeDb usable dynamic range (one-way equivalent × passes), dB
 * @param {number} o.attenuationDbPerKm fibre attenuation at the probe wavelength, dB/km
 * @param {number} [o.passes=2] @param {number} [o.fixedLossDb=0] one-way lumped losses, dB
 * @returns {number} sensing range, m
 */
export function sensingRange({ dynamicRangeDb, attenuationDbPerKm, passes = 2, fixedLossDb = 0 }) {
  pos(dynamicRangeDb, 'dynamicRangeDb'); pos(attenuationDbPerKm, 'attenuationDbPerKm'); pos(passes, 'passes'); nonNeg(fixedLossDb, 'fixedLossDb');
  return Math.max(0, ((dynamicRangeDb / passes - fixedLossDb) / attenuationDbPerKm) * 1000);
}

/**
 * Static-fatigue lifetime of a proof-tested fibre under a sustained service stress (power-law crack growth):
 * a fibre that survives the proof stress σp for time tp has a minimum time to failure
 *   ts = tp·(σp/σs)^n
 * under a constant stress σs (crack growth during proof-test unloading neglected), n being the stress-corrosion
 * susceptibility (about 20 for standard silica fibre; ≥ 18 required by fibre specifications).
 * Conversely the allowable sustained stress for a design life is σs = σp·(tp/life)^(1/n). Strain and stress are
 * related by σ = E·ε with E ≈ 72 GPa (first order); a 1 % proof strain is about 0.69 GPa (100 kpsi).
 * Reference: IEC TR 62048 (power-law reliability theory, simplified minimum-life form).
 * @param {object} o
 * @param {number} o.proofStress σp, Pa @param {number} [o.proofTime=1] tp, s
 * @param {number} [o.n=20] stress-corrosion parameter
 * @param {number} [o.serviceStress] σs, Pa @param {number} [o.designLife] s
 * @returns {{lifetime:?number, allowableStress:?number, allowableStressRatio:?number}} minimum lifetime at the
 *          service stress (s), allowable sustained stress for the design life (Pa) and its ratio to proof stress
 */
export function fibreStaticFatigue({ proofStress, proofTime = 1, n = 20, serviceStress, designLife }) {
  pos(proofStress, 'proofStress'); pos(proofTime, 'proofTime'); pos(n, 'n');
  let lifetime = null; let allowableStress = null;
  if (serviceStress !== undefined) { pos(serviceStress, 'serviceStress'); lifetime = proofTime * Math.pow(proofStress / serviceStress, n); }
  if (designLife !== undefined) { pos(designLife, 'designLife'); allowableStress = proofStress * Math.pow(proofTime / designLife, 1 / n); }
  return { lifetime, allowableStress, allowableStressRatio: allowableStress !== null ? allowableStress / proofStress : null };
}
