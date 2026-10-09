/**
 * Indicative property estimators for mooring-line components as functions of nominal diameter.
 *
 * These are preliminary-design estimators, not certified data: always replace them with the manufacturer's
 * certificate values for final design. Every coefficient is exposed (exported tables and per-call overrides).
 *
 * Sources:
 *  - Chain break load: IACS UR W22 / DNV-OS-E302 offshore mooring chain, BL = c d^2 (44 - 0.08 d) kN, d in mm.
 *  - Chain mass and axial stiffness, six-strand wire rope and 8-strand fibre rope coefficients: the widely used
 *    line-type estimators published in the OrcaFlex documentation ("Line type wizard": chain, rope/wire), which
 *    for wire rope coincide with the EN 12385-4 minimum breaking force factors at rope grade 1770.
 *  - Hydrodynamic coefficients and marine growth: DNV-OS-E301 Ch.2 Sec.2.
 *  - Fibre-rope stiffness model EA = (alpha + beta Tmean/MBL) MBL: Francois & Davies (2008); BV NI 432.
 */

const G = 9.80665;

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}

function posNum(v, name) {
  need(typeof v === 'number' && Number.isFinite(v) && v > 0, `${name} must be a finite number > 0`);
}

/** Break-load grade constants c in BL = c d^2 (44 - 0.08 d) [kN, d in mm] (IACS UR W22 / DNV-OS-E302). */
export const CHAIN_GRADE_C = Object.freeze({ R3: 0.0223, R3S: 0.0249, R4: 0.0274, R4S: 0.0304, R5: 0.0320 });

/**
 * Chain mass and stiffness coefficients (d in m): mass = massCoeff d^2 [kg/m], EA = eaCoeff d^2 [N].
 * Studless 19.9e3 kg/m per m^2 and 0.854e11 N per m^2; stud-link 21.9e3 and 1.01e11 (OrcaFlex chain estimators).
 */
export const CHAIN_COEFFS = Object.freeze({
  studless: Object.freeze({ massCoeff: 19.9e3, eaCoeff: 0.854e11 }),
  studlink: Object.freeze({ massCoeff: 21.9e3, eaCoeff: 1.01e11 }),
});

/**
 * Minimum breaking load of offshore mooring chain, BL = c d^2 (44 - 0.08 d) kN with d in mm (same for studless
 * and stud-link). Reference: IACS UR W22 Table (break load); DNV-OS-E302. Valid for the catalogued range of
 * about 50-210 mm; the formula peaks at d = 366 mm and must not be extrapolated towards it.
 * @param {object} o
 * @param {number} o.d nominal (bar) diameter [m]
 * @param {'R3'|'R3S'|'R4'|'R4S'|'R5'} [o.grade='R4']
 * @param {number} [o.c] override of the grade constant [-]
 * @returns {number} minimum breaking load [N]
 */
export function chainMBL({ d, grade = 'R4', c }) {
  posNum(d, 'd');
  need(d < 0.3, 'd must be below 0.3 m (formula range)');
  const cc = c ?? CHAIN_GRADE_C[grade];
  need(cc > 0, "grade must be one of 'R3', 'R3S', 'R4', 'R4S', 'R5' (or give c)");
  const dmm = d * 1000;
  return cc * dmm * dmm * (44 - 0.08 * dmm) * 1000;
}

function submerged(mass, rhoMat, rhoW, g) {
  const vol = mass / rhoMat; // displaced volume per metre [m^2]
  return { w: (mass - rhoW * vol) * g, dVol: Math.sqrt(4 * vol / Math.PI) };
}

/**
 * Indicative properties of studless or stud-link mooring chain.
 * mass = massCoeff d^2, EA = eaCoeff d^2 (CHAIN_COEFFS); MBL from chainMBL; submerged weight from the displaced
 * steel volume (first principles). Reference: see file header. Validity: preliminary design only.
 * @param {object} o
 * @param {number} o.d nominal diameter [m]
 * @param {'R3'|'R3S'|'R4'|'R4S'|'R5'} [o.grade='R4']
 * @param {'studless'|'studlink'} [o.type='studless']
 * @param {number} [o.massCoeff] override [kg/m per m^2] @param {number} [o.eaCoeff] override [N per m^2]
 * @param {number} [o.c] override of the break-load grade constant
 * @param {number} [o.rhoSteel=7850] steel density [kg/m^3] @param {number} [o.rhoWater=1025] [kg/m^3]
 * @param {number} [o.g=9.80665] [m/s^2]
 * @returns {{d:number,mass:number,w:number,EA:number,MBL:number,dVol:number,indicative:true}} mass in air
 *   [kg/m], submerged weight [N/m], EA [N], MBL [N], volume-equivalent diameter [m]
 */
export function chainProperties({ d, grade = 'R4', type = 'studless', massCoeff, eaCoeff, c, rhoSteel = 7850, rhoWater = 1025, g = G }) {
  posNum(d, 'd');
  need(type in CHAIN_COEFFS, "type must be 'studless' or 'studlink'");
  posNum(rhoSteel, 'rhoSteel'); posNum(rhoWater, 'rhoWater');
  const k = CHAIN_COEFFS[type];
  const mass = (massCoeff ?? k.massCoeff) * d * d;
  const s = submerged(mass, rhoSteel, rhoWater, g);
  return { d, mass, w: s.w, EA: (eaCoeff ?? k.eaCoeff) * d * d, MBL: chainMBL({ d, grade, c }), dVol: s.dVol, indicative: true };
}

/**
 * Wire-rope coefficients (d in m): mass = massCoeff d^2 [kg/m], EA = eaCoeff d^2 [N], MBL = mblCoeff d^2 [N].
 * Six-strand values are the OrcaFlex wire estimators (6x19 class, grade 1770; MBL = K d^2 Rr with the
 * EN 12385-4 factors K = 0.356 steel core and 0.330 fibre core). Spiral-strand values are indicative only:
 * typical ranges massCoeff 4.4e3-5.0e3, eaCoeff 8.5e10-10e10, mblCoeff 8.5e8-10e8.
 */
export const WIRE_COEFFS = Object.freeze({
  'six-strand-iwrc': Object.freeze({ massCoeff: 3.9897e3, eaCoeff: 4.04e10, mblCoeff: 6.33358e8 }),
  'six-strand-fc': Object.freeze({ massCoeff: 3.6109e3, eaCoeff: 3.67e10, mblCoeff: 5.84175e8 }),
  'spiral-strand': Object.freeze({ massCoeff: 4.69e3, eaCoeff: 9.0e10, mblCoeff: 9.0e8 }),
});

/**
 * Indicative properties of steel wire rope (six-strand with steel or fibre core, or spiral strand).
 * Submerged weight assumes the rope displaces its steel volume only (unsheathed, flooded voids); for sheathed
 * spiral strand pass dOuter so that buoyancy is taken on the sheathed envelope.
 * @param {object} o
 * @param {number} o.d nominal rope diameter [m]
 * @param {'six-strand-iwrc'|'six-strand-fc'|'spiral-strand'} [o.construction='six-strand-iwrc']
 * @param {number} [o.massCoeff] @param {number} [o.eaCoeff] @param {number} [o.mblCoeff] overrides (see WIRE_COEFFS)
 * @param {number} [o.dOuter] sheathed outer diameter [m]; buoyancy then uses the full envelope
 * @param {number} [o.rhoSteel=7850] [kg/m^3] @param {number} [o.rhoWater=1025] [kg/m^3] @param {number} [o.g]
 * @returns {{d:number,mass:number,w:number,EA:number,MBL:number,dVol:number,indicative:true}} see chainProperties
 */
export function wireRopeProperties({ d, construction = 'six-strand-iwrc', massCoeff, eaCoeff, mblCoeff, dOuter, rhoSteel = 7850, rhoWater = 1025, g = G }) {
  posNum(d, 'd');
  need(construction in WIRE_COEFFS, "construction must be 'six-strand-iwrc', 'six-strand-fc' or 'spiral-strand'");
  const k = WIRE_COEFFS[construction];
  const mass = (massCoeff ?? k.massCoeff) * d * d;
  let s;
  if (dOuter !== undefined) {
    need(dOuter >= d, 'dOuter must be >= d');
    s = { w: (mass - rhoWater * Math.PI * dOuter * dOuter / 4) * g, dVol: dOuter };
  } else s = submerged(mass, rhoSteel, rhoWater, g);
  return { d, mass, w: s.w, EA: (eaCoeff ?? k.eaCoeff) * d * d, MBL: (mblCoeff ?? k.mblCoeff) * d * d, dVol: s.dVol, indicative: true };
}

/**
 * Fibre-rope coefficients (d in m): mass = massCoeff d^2 [kg/m], EA = eaCoeff d^2 [N] (single linearised value),
 * MBL = mblCoeff d^2 [N], rhoFibre = fibre density [kg/m^3].
 * Polyester and nylon: 8-strand multiplait estimators from the OrcaFlex documentation (nylon MBL is the dry
 * value). HMPE: indicative mid-range values only (typical ranges massCoeff 0.5e3-0.65e3, mblCoeff 6e8-11e8,
 * eaCoeff 30-60 times MBL) - supply certificate values. Parallel-subrope deep-water polyester tethers are
 * considerably stronger per unit diameter than these 8-strand values.
 */
export const FIBRE_COEFFS = Object.freeze({
  polyester: Object.freeze({ massCoeff: 0.7978e3, eaCoeff: 1.09e9, mblCoeff: 1.70466e8, rhoFibre: 1380 }),
  nylon: Object.freeze({ massCoeff: 0.6476e3, eaCoeff: 1.18e8, mblCoeff: 1.39357e8, rhoFibre: 1140 }),
  hmpe: Object.freeze({ massCoeff: 0.58e3, eaCoeff: 3.8e10, mblCoeff: 8.5e8, rhoFibre: 975 }),
});

/**
 * Indicative properties of synthetic fibre rope. Submerged weight from the displaced fibre volume
 * w = m g (1 - rhoWater / rhoFibre) (first principles; HMPE floats, w < 0).
 * @param {object} o
 * @param {number} o.d nominal rope diameter [m]
 * @param {'polyester'|'nylon'|'hmpe'} [o.material='polyester']
 * @param {number} [o.massCoeff] @param {number} [o.eaCoeff] @param {number} [o.mblCoeff] @param {number} [o.rhoFibre]
 *   overrides (see FIBRE_COEFFS)
 * @param {number} [o.rhoWater=1025] [kg/m^3] @param {number} [o.g]
 * @returns {{d:number,mass:number,w:number,EA:number,MBL:number,dVol:number,indicative:true}} see chainProperties
 */
export function fibreRopeProperties({ d, material = 'polyester', massCoeff, eaCoeff, mblCoeff, rhoFibre, rhoWater = 1025, g = G }) {
  posNum(d, 'd');
  need(material in FIBRE_COEFFS, "material must be 'polyester', 'nylon' or 'hmpe'");
  const k = FIBRE_COEFFS[material];
  const mass = (massCoeff ?? k.massCoeff) * d * d;
  const s = submerged(mass, rhoFibre ?? k.rhoFibre, rhoWater, g);
  return { d, mass, w: s.w, EA: (eaCoeff ?? k.eaCoeff) * d * d, MBL: (mblCoeff ?? k.mblCoeff) * d * d, dVol: s.dVol, indicative: true };
}

/**
 * Presets for the fibre-rope stiffness model Kr = EA / MBL = alpha + beta (Tmean / MBL).
 * polyester-dynamic: Kr = 18.5 + 0.33 x mean load in % MBL (Francois & Davies 2008; BV NI 432), i.e.
 * alpha = 18.5, beta = 33 with the mean load as a fraction. polyester-static: a constant quasi-static value,
 * typical range 10-15 (rope and load-history dependent).
 */
export const FIBRE_STIFFNESS_PRESETS = Object.freeze({
  'polyester-dynamic': Object.freeze({ alpha: 18.5, beta: 33 }),
  'polyester-static': Object.freeze({ alpha: 13, beta: 0 }),
});

/**
 * Load-dependent axial stiffness of a fibre rope, EA = (alpha + beta Tmean / MBL) MBL.
 * Reference: Francois & Davies (2008) "Characterization of polyester mooring lines", OMAE2008-57136; BV NI 432.
 * Other materials or rope makes: supply alpha and beta from the rope data sheet (no default is assumed).
 * @param {object} o
 * @param {number} o.MBL minimum breaking load [N]
 * @param {number} o.Tmean mean tension [N] (0 <= Tmean <= MBL)
 * @param {number} [o.alpha] intercept [-] @param {number} [o.beta] slope per unit Tmean/MBL [-]
 * @param {'polyester-dynamic'|'polyester-static'} [o.preset] used when alpha / beta are not given
 * @returns {{EA:number,Kr:number,alpha:number,beta:number}} axial stiffness [N], Kr = EA/MBL [-]
 */
export function fibreRopeStiffness({ MBL, Tmean, alpha, beta, preset }) {
  posNum(MBL, 'MBL');
  need(Number.isFinite(Tmean) && Tmean >= 0 && Tmean <= MBL, 'Tmean must be within [0, MBL]');
  const p = preset ? FIBRE_STIFFNESS_PRESETS[preset] : undefined;
  need(!preset || p, "preset must be 'polyester-dynamic' or 'polyester-static'");
  const a = alpha ?? p?.alpha, b = beta ?? p?.beta;
  need(Number.isFinite(a) && a > 0, 'alpha must be given (> 0) or a preset selected');
  need(Number.isFinite(b) && b >= 0, 'beta must be given (>= 0) or a preset selected');
  const Kr = a + b * Tmean / MBL;
  return { EA: Kr * MBL, Kr, alpha: a, beta: b };
}

/**
 * Diameter after a corrosion / wear allowance: d_eff = d - rate x years. The break load of chain is then
 * re-evaluated on d_eff (chainMBL). DNV-OS-E301 Ch.2 Sec.2 tabulates allowances by zone (splash zone, catenary,
 * bottom), typically 0.2-0.8 mm/year; the rate is project/standard specific and must be supplied.
 * @param {object} o
 * @param {number} o.d new nominal diameter [m]
 * @param {number} o.rate diameter loss rate [m/year] (e.g. 0.4e-3)
 * @param {number} o.years exposure / design life [years]
 * @returns {number} reduced diameter [m]
 */
export function corrodedDiameter({ d, rate, years }) {
  posNum(d, 'd');
  need(Number.isFinite(rate) && rate >= 0, 'rate must be >= 0');
  need(Number.isFinite(years) && years >= 0, 'years must be >= 0');
  const de = d - rate * years;
  need(de > 0, 'rate x years exceeds the diameter');
  return de;
}

/**
 * Drag coefficients of mooring components, DNV-OS-E301 Ch.2 Sec.2 (table of transverse / longitudinal drag
 * coefficients; chain coefficients refer to the nominal bar diameter). Longitudinal values are tabulated for
 * chain only. Added-mass coefficients are not tabulated there: Can = 1.0 (per displaced volume) and
 * Cat = 0.5 chain / 0 rope are common analysis assumptions and are exposed as arguments.
 */
export const HYDRO_COEFFS = Object.freeze({
  'stud-chain': Object.freeze({ Cdn: 2.6, Cdt: 1.4 }),
  'studless-chain': Object.freeze({ Cdn: 2.4, Cdt: 1.15 }),
  'six-strand-rope': Object.freeze({ Cdn: 1.8, Cdt: null }),
  'spiral-rope-unsheathed': Object.freeze({ Cdn: 1.6, Cdt: null }),
  'spiral-rope-sheathed': Object.freeze({ Cdn: 1.2, Cdt: null }),
  'fibre-rope': Object.freeze({ Cdn: 1.6, Cdt: null }),
});

/**
 * Hydrodynamic coefficients for a mooring component.
 * Reference: DNV-OS-E301 Ch.2 Sec.2 (drag coefficients). Tangential drag of ropes and all added-mass
 * coefficients are analyst's choices: defaults are documented assumptions, not code values.
 * @param {object} o
 * @param {'stud-chain'|'studless-chain'|'six-strand-rope'|'spiral-rope-unsheathed'|'spiral-rope-sheathed'|'fibre-rope'} o.type
 * @param {number} [o.Cdn] @param {number} [o.Cdt] overrides [-]
 * @param {number} [o.Can=1.0] normal added-mass coefficient [-] (typical 1.0)
 * @param {number} [o.Cat] axial added-mass coefficient [-] (default 0.5 chain, 0 rope; typical 0-0.5)
 * @param {number} [o.ropeCdt=0] tangential drag coefficient used for ropes when Cdt is not given [-]
 *   (typical 0-0.05 on the pi d area)
 * @returns {{Cdn:number,Cdt:number,Can:number,Cat:number,source:string}}
 */
export function hydroCoefficients({ type, Cdn, Cdt, Can = 1.0, Cat, ropeCdt = 0 }) {
  const k = HYDRO_COEFFS[type];
  need(k, `type must be one of ${Object.keys(HYDRO_COEFFS).join(', ')}`);
  const chain = type.endsWith('chain');
  return {
    Cdn: Cdn ?? k.Cdn, Cdt: Cdt ?? k.Cdt ?? ropeCdt, Can, Cat: Cat ?? (chain ? 0.5 : 0),
    source: 'DNV-OS-E301 Ch.2 Sec.2 (drag); added mass and rope tangential drag are assumptions',
  };
}

/**
 * Marine growth on a mooring line, DNV-OS-E301 Ch.2 Sec.2 (marine growth):
 * mass  M = (pi/4) [(D + 2 t)^2 - D^2] rho_growth mu  [kg/m],
 * submerged weight W = M (1 - rho_sw / rho_growth) g  [N/m],
 * drag coefficient C_D,growth = C_D (D + 2 t) / D (still referred to the nominal diameter),
 * with mu = 2.0 for chain and 1.0 for wire / fibre rope.
 * @param {object} o
 * @param {number} o.d nominal diameter [m]
 * @param {number} o.thickness marine growth thickness t [m]
 * @param {number} [o.rhoGrowth=1325] marine growth density [kg/m^3] (site specific; 1325 is the NORSOK N-003 value)
 * @param {'chain'|'rope'} [o.component='chain'] sets mu (2.0 / 1.0)
 * @param {number} [o.mu] override [-]
 * @param {number} [o.Cdn] clean normal drag coefficient to be scaled [-]
 * @param {number} [o.rhoWater=1025] [kg/m^3] @param {number} [o.g]
 * @returns {{mass:number,w:number,CdnGrowth:number|null,dGrowth:number}} added mass in air [kg/m], added
 *   submerged weight [N/m], scaled drag coefficient [-], outer diameter with growth [m]
 */
export function marineGrowth({ d, thickness, rhoGrowth = 1325, component = 'chain', mu, Cdn, rhoWater = 1025, g = G }) {
  posNum(d, 'd');
  need(Number.isFinite(thickness) && thickness >= 0, 'thickness must be >= 0');
  posNum(rhoGrowth, 'rhoGrowth');
  need(component === 'chain' || component === 'rope', "component must be 'chain' or 'rope'");
  const m = mu ?? (component === 'chain' ? 2.0 : 1.0);
  const dg = d + 2 * thickness;
  const mass = (Math.PI / 4) * (dg * dg - d * d) * rhoGrowth * m;
  return { mass, w: mass * (1 - rhoWater / rhoGrowth) * g, CdnGrowth: Cdn !== undefined ? Cdn * dg / d : null, dGrowth: dg };
}

/**
 * Convenience: a multiseg.js segment record { L, w, EA, MBL, name } from an estimator result.
 * @param {object} props result of chainProperties / wireRopeProperties / fibreRopeProperties
 * @param {number} L unstretched length [m]
 * @param {string} [name]
 * @returns {{L:number,w:number,EA:number,MBL:number,mass:number,d:number,name:string}}
 */
export function toSegment(props, L, name = 'segment') {
  posNum(L, 'L');
  return { L, w: props.w, EA: props.EA, MBL: props.MBL, mass: props.mass, d: props.dVol, name };
}
