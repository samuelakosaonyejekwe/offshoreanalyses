// Armour stability for rubble-mound structures: Hudson, Van der Meer, concrete units, toe, gradings,
// filters and cross-section quantities. SI units (masses in kg, densities in kg/m³).

const G = 9.80665;

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function pos(v, name) {
  need(Number.isFinite(v) && v > 0, `${name} must be a positive finite number`);
}

/**
 * Relative buoyant density Δ = ρs/ρw − 1.
 * @param {object} o
 * @param {number} o.rhoS armour density [kg/m³]
 * @param {number} [o.rhoW] water density [kg/m³], default 1025
 * @returns {number} Δ [-]
 */
export function relativeDensity({ rhoS, rhoW = 1025 }) {
  pos(rhoS, 'rhoS');
  pos(rhoW, 'rhoW');
  need(rhoS > rhoW, 'rhoS must exceed rhoW');
  return rhoS / rhoW - 1;
}

/**
 * Hudson stability coefficients KD for no-damage design (0–5 % displacement), structure trunk.
 * Each entry gives the value for breaking and non-breaking waves at the structure.
 * Reference: Shore Protection Manual (1984) Table 7-8 (n = 2 layers, random placement unless noted).
 * Units: [-].
 */
export const HUDSON_KD = Object.freeze({
  'quarrystone-smooth-rounded': { breaking: 1.2, nonBreaking: 2.4 },
  'quarrystone-rough-angular': { breaking: 2.0, nonBreaking: 4.0 },
  'quarrystone-rough-angular-3plus-layers': { breaking: 2.2, nonBreaking: 4.5 },
  tetrapod: { breaking: 7.0, nonBreaking: 8.0 },
  quadripod: { breaking: 7.0, nonBreaking: 8.0 },
  tribar: { breaking: 9.0, nonBreaking: 10.0 },
  dolos: { breaking: 15.8, nonBreaking: 31.8 },
  'modified-cube': { breaking: 6.5, nonBreaking: 7.5 },
});

/**
 * Hudson (1959) armour mass: M50 = ρs H³ / (KD Δ³ cot α); equivalently Ns = H/(Δ Dn50) = (KD cot α)^(1/3).
 * Reference: Hudson (1959) J. Waterways Harbors Div. ASCE 85(WW3); SPM (1984) Eq. 7-116
 * (SPM uses H = H1/10 for random waves; CEM/Rock Manual often use Hs with adjusted KD).
 * Validity: 1.5 <= cot α <= 4 (SPM), non-overtopped permeable-core structures, no period or duration effect.
 * @param {object} o
 * @param {number} o.H design wave height at the structure [m]
 * @param {number} o.KD stability coefficient [-] (see HUDSON_KD)
 * @param {number} o.cotAlpha slope cotangent [-]
 * @param {number} o.rhoS armour density [kg/m³]
 * @param {number} [o.rhoW] water density [kg/m³]
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{M50:number, W50:number, Dn50:number, Ns:number, delta:number}} mass [kg], weight [N],
 *   nominal diameter [m], stability number [-], Δ [-]
 */
export function hudson({ H, KD, cotAlpha, rhoS, rhoW = 1025, g = G }) {
  pos(H, 'H');
  pos(KD, 'KD');
  pos(cotAlpha, 'cotAlpha');
  const delta = relativeDensity({ rhoS, rhoW });
  const M50 = (rhoS * H * H * H) / (KD * delta * delta * delta * cotAlpha);
  const Dn50 = Math.cbrt(M50 / rhoS);
  return { M50, W50: M50 * g, Dn50, Ns: Math.cbrt(KD * cotAlpha), delta };
}

function vdmCore({ xi, cotAlpha, P, sn, cpl, cs }) {
  const tanA = 1 / cotAlpha;
  const xiCr = Math.pow((cpl / cs) * Math.pow(P, 0.31) * Math.sqrt(tanA), 1 / (P + 0.5));
  const plunging = cpl * Math.pow(P, 0.18) * Math.pow(sn, 0.2) / Math.sqrt(xi);
  const surging = cs * Math.pow(P, -0.13) * Math.pow(sn, 0.2) * Math.sqrt(cotAlpha) * Math.pow(xi, P);
  // For cot α >= 4 the transition to surging waves does not exist and only the plunging formula applies.
  const usePl = xi < xiCr || cotAlpha >= 4;
  return { xiCr, plunging, surging, Ns: usePl ? plunging : surging, regime: usePl ? 'plunging' : 'surging' };
}

function vdmInputs({ Hs, T, cotAlpha, P, S, N, g }) {
  pos(Hs, 'Hs');
  pos(T, 'T');
  pos(cotAlpha, 'cotAlpha');
  need(Number.isFinite(P) && P >= 0.1 && P <= 0.6, 'P must lie in [0.1, 0.6]');
  pos(S, 'S');
  pos(N, 'N');
  const L0 = (g * T * T) / (2 * Math.PI);
  return { xi: 1 / cotAlpha / Math.sqrt(Hs / L0), sn: S / Math.sqrt(Math.min(N, 7500)) };
}

/**
 * Van der Meer (1988) stability of rock armour in deep water:
 *   plunging (ξm < ξcr): Hs/(Δ Dn50) = 6.2 P^0.18 (S/√N)^0.2 ξm^−0.5
 *   surging  (ξm >= ξcr): Hs/(Δ Dn50) = 1.0 P^−0.13 (S/√N)^0.2 √cot α ξm^P
 *   ξcr = [6.2 P^0.31 √tan α]^(1/(P+0.5));  for cot α >= 4 only the plunging formula is used.
 * Reference: Van der Meer (1988) "Rock slopes and gravel beaches under wave attack", Delft Hydraulics
 * Publ. 396; Rock Manual (CIRIA C683, 2007) §5.2.2.2.
 * Validity: 0.1 <= P <= 0.6, 0.005 <= sm <= 0.06, N <= 7500 (damage is capped at N = 7500), 1.5 <= cot α <= 6,
 * deep water at the toe (h > 3 Hs).
 * @param {object} o
 * @param {number} o.Hs significant wave height at the toe [m]
 * @param {number} o.Tm mean wave period [s]
 * @param {number} o.cotAlpha slope cotangent [-]
 * @param {number} o.P notional permeability [-] (0.1 impermeable core … 0.6 homogeneous)
 * @param {number} o.S damage level [-] (2–3 start of damage, 8–17 failure depending on slope)
 * @param {number} o.N number of waves [-]
 * @param {number} o.delta relative buoyant density Δ [-]
 * @param {number} [o.rhoS] armour density for the mass [kg/m³]
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{Ns:number,Dn50:number,M50:number,xi:number,xiCr:number,regime:string}} stability number [-],
 *   Dn50 [m], M50 [kg] (NaN without rhoS), ξm [-], ξcr [-], 'plunging' | 'surging'
 */
export function vanDerMeerRock({ Hs, Tm, cotAlpha, P, S, N, delta, rhoS, g = G }) {
  pos(delta, 'delta');
  const { xi, sn } = vdmInputs({ Hs, T: Tm, cotAlpha, P, S, N, g });
  const r = vdmCore({ xi, cotAlpha, P, sn, cpl: 6.2, cs: 1.0 });
  const Dn50 = Hs / (delta * r.Ns);
  return { Ns: r.Ns, Dn50, M50: rhoS ? rhoS * Dn50 ** 3 : Number.NaN, xi, xiCr: r.xiCr, regime: r.regime };
}

/**
 * Damage level S of a given rock armour from the Van der Meer (1988) formulae (inverse of vanDerMeerRock).
 * Reference and validity: as vanDerMeerRock.
 * @param {object} o
 * @param {number} o.Hs [m]
 * @param {number} o.Tm [s]
 * @param {number} o.cotAlpha [-]
 * @param {number} o.P [-]
 * @param {number} o.N number of waves [-]
 * @param {number} o.delta Δ [-]
 * @param {number} o.Dn50 nominal diameter [m]
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{S:number, Ns:number, regime:string}} damage level [-], Hs/(Δ Dn50) [-], regime
 */
export function vanDerMeerDamage({ Hs, Tm, cotAlpha, P, N, delta, Dn50, g = G }) {
  pos(delta, 'delta');
  pos(Dn50, 'Dn50');
  const { xi } = vdmInputs({ Hs, T: Tm, cotAlpha, P, S: 1, N, g });
  // With S/√N = 1 the core returns the bracket multiplying (S/√N)^0.2.
  const r = vdmCore({ xi, cotAlpha, P, sn: 1, cpl: 6.2, cs: 1.0 });
  const Ns = Hs / (delta * Dn50);
  return { S: Math.sqrt(Math.min(N, 7500)) * Math.pow(Ns / r.Ns, 5), Ns, regime: r.regime };
}

/**
 * Van der Meer formulae modified for shallow water (Van Gent et al. 2003):
 *   plunging: Hs/(Δ Dn50) = cpl P^0.18 (S/√N)^0.2 (Hs/H2%) ξs−1,0^−0.5
 *   surging:  Hs/(Δ Dn50) = cs P^−0.13 (S/√N)^0.2 (Hs/H2%) √cot α ξs−1,0^P
 *   cpl = 8.4, cs = 1.3, ξcr = [(cpl/cs) P^0.31 √tan α]^(1/(P+0.5)), ξs−1,0 from Hs and Tm−1,0 at the toe.
 * Reference: Van Gent, Smale & Kuiper (2003) Proc. Coastal Structures 2003; Rock Manual (2007)
 * §5.2.2.2.
 * Validity: shallow foreshores (about 1.5 < h/Hs,toe < 3), N < 3000; coefficients exposed as arguments.
 * @param {object} o
 * @param {number} o.Hs significant wave height H1/3 at the toe [m]
 * @param {number} o.H2 wave height exceeded by 2 % of the waves at the toe [m]
 * @param {number} o.Tm10 spectral period Tm−1,0 at the toe [s]
 * @param {number} o.cotAlpha [-]
 * @param {number} o.P [-]
 * @param {number} o.S [-]
 * @param {number} o.N [-]
 * @param {number} o.delta [-]
 * @param {number} [o.cpl] default 8.4
 * @param {number} [o.cs] default 1.3
 * @param {number} [o.rhoS] [kg/m³]
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{Ns:number,Dn50:number,M50:number,xi:number,xiCr:number,regime:string}} as vanDerMeerRock
 */
export function vanDerMeerShallow({ Hs, H2, Tm10, cotAlpha, P, S, N, delta, cpl = 8.4, cs = 1.3, rhoS, g = G }) {
  pos(delta, 'delta');
  pos(H2, 'H2');
  const { xi, sn } = vdmInputs({ Hs, T: Tm10, cotAlpha, P, S, N, g });
  const r = vdmCore({ xi, cotAlpha, P, sn, cpl, cs });
  const Ns = r.Ns * (Hs / H2);
  const Dn50 = Hs / (delta * Ns);
  return { Ns, Dn50, M50: rhoS ? rhoS * Dn50 ** 3 : Number.NaN, xi, xiCr: r.xiCr, regime: r.regime };
}

/**
 * Van der Meer (1988) stability of randomly placed concrete units in a double layer on a 1:1.5 slope:
 *   cubes:     Hs/(Δ Dn) = (6.7 Nod^0.4 / N^0.3 + 1.0) som^−0.1
 *   tetrapods: Hs/(Δ Dn) = (3.75 Nod^0.5 / N^0.25 + 0.85) som^−0.2
 * with som = 2π Hs/(g Tm²) and Nod the number of displaced units in a strip one Dn wide.
 * Reference: Van der Meer (1988) Proc. 21st ICCE "Stability of cubes, tetrapods and Accropode";
 * Rock Manual (2007) §5.2.2.3.
 * Validity: slope 1:1.5, non-depth-limited waves, surging conditions (3 < ξm < 6); Nod = 0 gives
 * the no-damage criterion.
 * @param {object} o
 * @param {'cube'|'tetrapod'} o.unit
 * @param {number} o.Hs [m]
 * @param {number} o.Tm mean period [s]
 * @param {number} o.Nod damage number [-] (0 no damage; 0.5 start; 2 failure for cubes, 1.5 for tetrapods)
 * @param {number} o.N number of waves [-]
 * @param {number} o.delta Δ [-]
 * @param {number} [o.rhoS] concrete density [kg/m³]
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{Ns:number, Dn:number, M:number, som:number}} stability number [-], nominal diameter [m],
 *   unit mass [kg] (NaN without rhoS), steepness [-]
 */
export function concreteUnitVanDerMeer({ unit, Hs, Tm, Nod, N, delta, rhoS, g = G }) {
  pos(Hs, 'Hs');
  pos(Tm, 'Tm');
  pos(N, 'N');
  pos(delta, 'delta');
  need(Number.isFinite(Nod) && Nod >= 0, 'Nod must be >= 0');
  const som = (2 * Math.PI * Hs) / (g * Tm * Tm);
  let Ns;
  if (unit === 'cube') Ns = ((6.7 * Math.pow(Nod, 0.4)) / Math.pow(N, 0.3) + 1.0) * Math.pow(som, -0.1);
  else if (unit === 'tetrapod') Ns = ((3.75 * Math.sqrt(Nod)) / Math.pow(N, 0.25) + 0.85) * Math.pow(som, -0.2);
  else throw new RangeError("unit must be 'cube' or 'tetrapod'");
  const Dn = Hs / (delta * Ns);
  return { Ns, Dn, M: rhoS ? rhoS * Dn ** 3 : Number.NaN, som };
}

/**
 * Unit size from a design stability number Ns = Hs/(Δ Dn) (e.g. single-layer interlocking units, for which
 * the licensor's design value applies; the number is an argument because it is unit- and slope-specific).
 * Reference: definition of the stability number, Rock Manual (2007) §5.2.2.3.
 * @param {object} o
 * @param {number} o.Hs [m]
 * @param {number} o.Ns design stability number [-] (typical published design values 2.5–2.8 for
 *   single-layer interlocking units)
 * @param {number} o.delta Δ [-]
 * @param {number} [o.rhoS] unit density [kg/m³]
 * @returns {{Dn:number, M:number, volume:number}} nominal diameter [m], mass [kg], unit volume [m³]
 */
export function stabilityNumberSize({ Hs, Ns, delta, rhoS }) {
  pos(Hs, 'Hs');
  pos(Ns, 'Ns');
  pos(delta, 'delta');
  const Dn = Hs / (delta * Ns);
  return { Dn, M: rhoS ? rhoS * Dn ** 3 : Number.NaN, volume: Dn ** 3 };
}

/**
 * Toe-berm rock stability, Van der Meer (1998): Hs/(Δ Dn50) = (2 + 6.2 (ht/h)^2.7) Nod^0.15.
 * Reference: Van der Meer (1998) in "Dikes and revetments" (ed. Pilarczyk), ch. 11;
 * Rock Manual (2007) §5.2.2.9.
 * Validity: 0.4 < ht/h < 0.9, 3 < ht/Dn50 < 25, depth-limited conditions; Nod = 0.5 start of damage,
 * 2 some flattening, 4 complete flattening of the toe.
 * @param {object} o
 * @param {number} o.Hs significant wave height at the toe [m]
 * @param {number} o.ht water depth above the toe berm [m]
 * @param {number} o.h water depth in front of the toe [m]
 * @param {number} o.delta Δ [-]
 * @param {number} [o.Nod] damage number [-], default 2
 * @param {number} [o.rhoS] [kg/m³]
 * @returns {{Ns:number,Dn50:number,M50:number,valid:boolean}} stability number [-], Dn50 [m], M50 [kg],
 *   validity of ht/h and ht/Dn50
 */
export function toeStabilityVanDerMeer({ Hs, ht, h, delta, Nod = 2, rhoS }) {
  pos(Hs, 'Hs');
  pos(ht, 'ht');
  pos(h, 'h');
  pos(delta, 'delta');
  pos(Nod, 'Nod');
  need(ht <= h, 'ht must not exceed h');
  const Ns = (2 + 6.2 * Math.pow(ht / h, 2.7)) * Math.pow(Nod, 0.15);
  const Dn50 = Hs / (delta * Ns);
  const r = ht / h;
  const valid = r > 0.4 && r < 0.9 && ht / Dn50 > 3 && ht / Dn50 < 25;
  return { Ns, Dn50, M50: rhoS ? rhoS * Dn50 ** 3 : Number.NaN, valid };
}

/**
 * Nominal diameter from mass: Dn50 = (M50/ρs)^(1/3); sieve diameter D50 = Dn50/shape (shape ≈ 0.84).
 * Reference: Rock Manual (2007) §3.4.2 (Dn50 ≈ 0.84 D50).
 * @param {object} o
 * @param {number} o.M50 median mass [kg]
 * @param {number} o.rhoS rock density [kg/m³]
 * @param {number} [o.shape] Dn50/D50 [-], default 0.84
 * @returns {{Dn50:number, D50:number}} [m]
 */
export function massToDn({ M50, rhoS, shape = 0.84 }) {
  pos(M50, 'M50');
  pos(rhoS, 'rhoS');
  pos(shape, 'shape');
  const Dn50 = Math.cbrt(M50 / rhoS);
  return { Dn50, D50: Dn50 / shape };
}

/**
 * Median mass from nominal diameter: M50 = ρs Dn50³.
 * Reference: Rock Manual (2007) §3.4.2.
 * @param {object} o
 * @param {number} o.Dn50 [m]
 * @param {number} o.rhoS [kg/m³]
 * @returns {number} M50 [kg]
 */
export function dnToMass({ Dn50, rhoS }) {
  pos(Dn50, 'Dn50');
  pos(rhoS, 'rhoS');
  return rhoS * Dn50 ** 3;
}

/**
 * Armour layer thickness t = n kt Dn50, placing density N = n kt (1 − nv)/Dn50² and placed mass per area.
 * Reference: SPM (1984) Eq. 7-121, 7-122; Rock Manual (2007) §3.5.1.
 * Validity: kt (layer coefficient) and nv (layer porosity) depend on unit type and placement —
 * typical rough quarrystone, 2 layers random: kt = 1.0, nv = 0.37.
 * @param {object} o
 * @param {number} o.Dn50 [m]
 * @param {number} [o.n] number of layers [-], default 2
 * @param {number} [o.kt] layer coefficient [-], default 1.0
 * @param {number} [o.porosity] layer porosity nv [-], default 0.37
 * @param {number} [o.rhoS] density for the placed mass [kg/m³]
 * @returns {{thickness:number, unitsPerArea:number, massPerArea:number}} [m], [1/m²], [kg/m²]
 */
export function armourLayer({ Dn50, n = 2, kt = 1.0, porosity = 0.37, rhoS }) {
  pos(Dn50, 'Dn50');
  pos(n, 'n');
  pos(kt, 'kt');
  need(porosity >= 0 && porosity < 1, 'porosity must lie in [0, 1)');
  const thickness = n * kt * Dn50;
  const unitsPerArea = (n * kt * (1 - porosity)) / (Dn50 * Dn50);
  return { thickness, unitsPerArea, massPerArea: rhoS ? unitsPerArea * rhoS * Dn50 ** 3 : Number.NaN };
}

/**
 * Geometrically closed granular filter check between a filter (upper, coarser) layer and its base:
 *   retention     D15f/D85b <= retentionLimit        (Terzaghi; 4–5)
 *   permeability  D15f/D15b >= permeabilityLimit     (Terzaghi; 4–5)
 *   uniformity    D50f/D50b <= uniformityLimit       (USACE; 25)
 * and the breakwater underlayer rule M50,under >= M50,armour/massRatio (SPM 1984: 10; Rock Manual: 10–15).
 * Reference: Terzaghi & Peck (1948); SPM (1984) §7.III.7.g; Rock Manual (2007) §5.4.3.6.
 * Validity: stable (non-gap-graded) gradings under cyclic flow; limits are arguments.
 * @param {object} o
 * @param {number} o.D15f [m]
 * @param {number} [o.D50f] [m]
 * @param {number} o.D85b [m]
 * @param {number} o.D15b [m]
 * @param {number} [o.D50b] [m]
 * @param {number} [o.retentionLimit] default 5
 * @param {number} [o.permeabilityLimit] default 5
 * @param {number} [o.uniformityLimit] default 25
 * @returns {{retention:{ratio:number,ok:boolean},permeability:{ratio:number,ok:boolean},
 *   uniformity:{ratio:number,ok:boolean},ok:boolean}} ratios [-]
 */
export function filterRules({ D15f, D50f, D85b, D15b, D50b, retentionLimit = 5, permeabilityLimit = 5, uniformityLimit = 25 }) {
  pos(D15f, 'D15f');
  pos(D85b, 'D85b');
  pos(D15b, 'D15b');
  const r = D15f / D85b;
  const p = D15f / D15b;
  const hasU = D50f !== undefined && D50b !== undefined;
  const u = hasU ? D50f / D50b : Number.NaN;
  const retention = { ratio: r, ok: r <= retentionLimit };
  const permeability = { ratio: p, ok: p >= permeabilityLimit };
  const uniformity = { ratio: u, ok: hasU ? u <= uniformityLimit : true };
  return { retention, permeability, uniformity, ok: retention.ok && permeability.ok && uniformity.ok };
}

/**
 * Underlayer mass from the armour mass: M50,under = M50,armour/ratio.
 * Reference: SPM (1984) §7.III.7 (W/10); Rock Manual (2007) §5.2.2.2 (1/10 to 1/15).
 * @param {object} o
 * @param {number} o.M50armour [kg]
 * @param {number} [o.ratio] default 10 (range 10–15)
 * @returns {number} M50 of the underlayer [kg]
 */
export function underlayerMass({ M50armour, ratio = 10 }) {
  pos(M50armour, 'M50armour');
  pos(ratio, 'ratio');
  return M50armour / ratio;
}

/**
 * Quantity take-off for a trapezoidal rubble-mound cross-section built of concentric layers.
 * The outer profile is a trapezoid of crest width B, height Hh above the seabed and slopes cot α_front,
 * cot α_rear. Each layer has thickness t measured normal to the slopes and vertically on the crest, so the
 * next profile inside is a trapezoid of height Hh − t and crest width
 *   B − t [(√(1+cf²) − cf) + (√(1+cr²) − cr)]      (offset geometry, first principles)
 * and the area of a trapezoid is A = B Hh + ½ Hh² (cf + cr). What remains inside the last layer is the core.
 * Reference: plane geometry; bulk-to-placed conversion with layer porosity per Rock Manual (2007) §3.5.
 * Validity: flat seabed, layers wrap front slope, crest and rear slope with constant thickness.
 * @param {object} o
 * @param {number} o.crestWidth B [m]
 * @param {number} o.height Hh, crest level above seabed [m]
 * @param {number} o.cotFront front slope cotangent [-]
 * @param {number} o.cotRear rear slope cotangent [-]
 * @param {Array<{name:string,thickness:number,porosity?:number,rhoS?:number}>} o.layers outer to inner;
 *   thickness [m], porosity [-], rhoS [kg/m³]
 * @param {{porosity?:number, rhoS?:number}} [o.core] core porosity [-] and density [kg/m³]
 * @returns {{totalArea:number, baseWidth:number, layers:Array<{name:string,area:number,solidVolume:number,
 *   mass:number,crestWidth:number,height:number}>}} areas = bulk volume per metre run [m³/m], solid volume
 *   [m³/m], mass [kg/m] (NaN without rhoS); last entry is the core; crestWidth/height are the outer
 *   dimensions of that layer [m]
 */
export function rubbleMoundQuantities({ crestWidth, height, cotFront, cotRear, layers = [], core = {} }) {
  pos(crestWidth, 'crestWidth');
  pos(height, 'height');
  need(Number.isFinite(cotFront) && cotFront >= 0, 'cotFront must be >= 0');
  need(Number.isFinite(cotRear) && cotRear >= 0, 'cotRear must be >= 0');
  const cs = cotFront + cotRear;
  const shrink = Math.sqrt(1 + cotFront * cotFront) - cotFront + Math.sqrt(1 + cotRear * cotRear) - cotRear;
  const area = (B, H) => B * H + 0.5 * H * H * cs;
  let B = crestWidth;
  let H = height;
  const totalArea = area(B, H);
  const out = [];
  const row = (name, a, por = 0, rhoS, Bo, Ho) => ({
    name, area: a, solidVolume: a * (1 - por), mass: rhoS ? a * (1 - por) * rhoS : Number.NaN, crestWidth: Bo, height: Ho,
  });
  for (const l of layers) {
    pos(l.thickness, `layers[${l.name}].thickness`);
    const Bi = B - l.thickness * shrink;
    const Hi = H - l.thickness;
    need(Bi > 0 && Hi > 0, `layer ${l.name} does not fit inside the cross-section`);
    out.push(row(l.name, area(B, H) - area(Bi, Hi), l.porosity, l.rhoS, B, H));
    B = Bi;
    H = Hi;
  }
  out.push(row('core', area(B, H), core.porosity, core.rhoS, B, H));
  return { totalArea, baseWidth: crestWidth + height * cs, layers: out };
}
