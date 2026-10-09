// Layered soil profile for pile design: property lookup, vertical effective stress, scour.
// Depth z is measured downwards from the ORIGINAL (unscoured) mudline, in metres. All values SI.

const G = 9.80665;
const TYPES = ['clay', 'sand', 'silt', 'rock'];

function num(v) {
  return typeof v === 'number' && Number.isFinite(v);
}

/**
 * Check a list of layers without throwing.
 * Layer fields (all SI): type 'clay'|'sand'|'silt'|'rock'; top [m]; thickness [m]; gammaSub [N/m³] submerged unit
 * weight; suTop [Pa] and suGrad [Pa/m] undrained shear strength at layer top and its gradient; phiDeg [deg]
 * effective friction angle; Dr [-] relative density 0..1; eps50 [-]; OCR [-]; G0 [Pa] small-strain shear modulus
 * (or Vs [m/s] shear-wave velocity); k [N/m³] initial modulus of subgrade reaction; PI [%] plasticity index;
 * St [-] sensitivity. Any further fields (e.g. apiSand, qur, Erock, J) are carried through unchanged.
 * Reference: none (data model). Validity: layers must be contiguous from z = 0 without gaps or overlaps.
 * @param {object[]} layers layer records as described above
 * @returns {{ valid: boolean, errors: string[] }} validation outcome with one message per defect
 */
export function validateLayers(layers) {
  const errors = [];
  if (!Array.isArray(layers) || layers.length === 0) return { valid: false, errors: ['layers: must be a non-empty array'] };
  const sorted = layers.slice().sort((a, b) => (a.top ?? 0) - (b.top ?? 0));
  let expectTop = 0;
  sorted.forEach((l, i) => {
    const tag = `layers[${i}]`;
    if (!TYPES.includes(l.type)) errors.push(`${tag}.type: must be one of ${TYPES.join(', ')}`);
    if (!num(l.top) || l.top < 0) errors.push(`${tag}.top: must be a number >= 0`);
    if (!num(l.thickness) || l.thickness <= 0) errors.push(`${tag}.thickness: must be > 0`);
    if (!num(l.gammaSub) || l.gammaSub <= 0 || l.gammaSub > 30e3) errors.push(`${tag}.gammaSub: must be in (0, 30e3] N/m³`);
    if (num(l.top) && Math.abs(l.top - expectTop) > 1e-9) errors.push(`${tag}.top: expected ${expectTop} (layers must be contiguous from 0)`);
    if (num(l.top) && num(l.thickness)) expectTop = l.top + l.thickness;
    if (l.suTop !== undefined && (!num(l.suTop) || l.suTop < 0)) errors.push(`${tag}.suTop: must be >= 0 Pa`);
    if (l.suGrad !== undefined && !num(l.suGrad)) errors.push(`${tag}.suGrad: must be a number`);
    if (num(l.suTop) && num(l.suGrad) && num(l.thickness) && l.suTop + l.suGrad * l.thickness < 0) errors.push(`${tag}.suGrad: su becomes negative within the layer`);
    if (l.phiDeg !== undefined && (!num(l.phiDeg) || l.phiDeg <= 0 || l.phiDeg >= 50)) errors.push(`${tag}.phiDeg: must be in (0, 50) deg`);
    if (l.Dr !== undefined && (!num(l.Dr) || l.Dr < 0 || l.Dr > 1)) errors.push(`${tag}.Dr: must be in [0, 1]`);
    if (l.eps50 !== undefined && (!num(l.eps50) || l.eps50 <= 0 || l.eps50 > 0.1)) errors.push(`${tag}.eps50: must be in (0, 0.1]`);
    if (l.OCR !== undefined && (!num(l.OCR) || l.OCR < 1)) errors.push(`${tag}.OCR: must be >= 1`);
    for (const key of ['G0', 'Vs', 'k', 'St']) {
      if (l[key] !== undefined && (!num(l[key]) || l[key] <= 0)) errors.push(`${tag}.${key}: must be > 0`);
    }
    if (l.PI !== undefined && (!num(l.PI) || l.PI < 0)) errors.push(`${tag}.PI: must be >= 0 %`);
    if (l.type === 'clay' && l.suTop === undefined) errors.push(`${tag}.suTop: required for clay`);
    if (l.type === 'sand' && l.phiDeg === undefined) errors.push(`${tag}.phiDeg: required for sand`);
  });
  return { valid: errors.length === 0, errors };
}

/**
 * Build a normalised soil profile from layer records (see validateLayers for the fields).
 * Reference: none (data model). Validity: contiguous layers from the mudline.
 * @param {object[]} layers layer records, SI units
 * @param {{ gammaW?: number }} [opts] gammaW [N/m³] unit weight of pore water (default 10 052, sea water)
 * @returns {{ layers: object[], depth: number, gammaW: number, scour: null }} profile; depth [m] is the base of
 *   the last layer; each layer gains bottom [m] and sigmaTop [Pa] (unscoured vertical effective stress at its top)
 */
export function createProfile(layers, opts = {}) {
  const v = validateLayers(layers);
  if (!v.valid) throw new RangeError(`layers: ${v.errors.join('; ')}`);
  const gammaW = opts.gammaW ?? 10052;
  if (!num(gammaW) || gammaW <= 0) throw new RangeError('gammaW: must be > 0 N/m³');
  let sigma = 0;
  const out = layers.slice().sort((a, b) => a.top - b.top).map((l) => {
    const layer = { suGrad: 0, ...l, bottom: l.top + l.thickness, sigmaTop: sigma };
    sigma += l.gammaSub * l.thickness;
    return layer;
  });
  return { layers: out, depth: out[out.length - 1].bottom, gammaW, scour: null };
}

/**
 * Index of the layer containing depth z. A depth exactly on an interface belongs to the lower layer.
 * @param {object} profile profile from createProfile
 * @param {number} z depth below original mudline [m], 0 <= z <= profile.depth
 * @returns {number} layer index [-]
 */
export function layerIndexAt(profile, z) {
  if (!num(z) || z < -1e-9 || z > profile.depth + 1e-9) throw new RangeError(`z: ${z} m is outside the profile (0..${profile.depth} m)`);
  const ls = profile.layers;
  for (let i = 0; i < ls.length; i++) if (z < ls[i].bottom - 1e-12) return i;
  return ls.length - 1;
}

function sigmaUnscoured(profile, z) {
  const l = profile.layers[layerIndexAt(profile, z)];
  return l.sigmaTop + l.gammaSub * (Math.max(z, 0) - l.top);
}

/**
 * Depth of the soil surface after scour (base of the local scour hole), measured from the original mudline.
 * @param {object} profile profile, possibly scoured
 * @returns {number} depth of the scoured surface at the pile [m] (0 when no scour)
 */
export function scourBaseDepth(profile) {
  return profile.scour ? profile.scour.general + profile.scour.local : 0;
}

/**
 * Vertical effective stress σ′v(z) = ∫ γ′ dz, including the scour adjustment if the profile carries one.
 * With scour: σ′v = 0 above the local scour base; below the overburden reduction depth it equals the
 * general-scour value σ′v,unscoured(z) − σ′v,unscoured(general); in between it varies linearly from zero at the
 * local scour base to the general-scour value at the reduction depth.
 * Reference: API RP 2A-WSD 21st ed. commentary on scour (C6.3) as carried into API RP 2GEO (2011); ISO 19901-4.
 * Validity: hydrostatic pore pressure, 0 <= z <= profile.depth.
 * @param {object} profile profile from createProfile / applyScour
 * @param {number} z depth below original mudline [m]
 * @returns {number} vertical effective stress [Pa]
 */
export function verticalEffectiveStress(profile, z) {
  const s0 = sigmaUnscoured(profile, z);
  const sc = profile.scour;
  if (!sc) return s0;
  const zl = sc.general + sc.local;
  if (z <= zl) return 0;
  const sGen = (zz) => sigmaUnscoured(profile, zz) - sigmaUnscoured(profile, sc.general);
  const zr = sc.general + sc.reductionDepth;
  if (z >= zr || zr <= zl) return sGen(z);
  return sGen(zr) * (z - zl) / (zr - zl);
}

/**
 * Soil properties at a depth. Linear su within a layer; G0 from Vs as ρ·Vs² when G0 is not given, with the
 * saturated density ρ = (γ′ + γw)/g.
 * Reference: none (lookup); G0 = ρ Vs² is the elastic shear-wave relation.
 * @param {object} profile profile from createProfile / applyScour
 * @param {number} z depth below original mudline [m]
 * @returns {object} { layerIndex, type, z [m], zBelowSurface [m] (below scoured surface, >= 0), sigmaV [Pa],
 *   gammaSub [N/m³], su [Pa]|undefined, phiDeg, Dr, eps50, OCR, G0 [Pa]|undefined, Vs [m/s]|undefined,
 *   k [N/m³]|undefined, PI [%], St, layer (the full layer record incl. pass-through fields) }
 */
export function propertiesAt(profile, z) {
  const i = layerIndexAt(profile, z);
  const l = profile.layers[i];
  const su = l.suTop === undefined ? undefined : l.suTop + l.suGrad * (z - l.top);
  let G0 = l.G0;
  if (G0 === undefined && l.Vs !== undefined) G0 = ((l.gammaSub + profile.gammaW) / G) * l.Vs * l.Vs;
  return {
    layerIndex: i, type: l.type, z, zBelowSurface: Math.max(z - scourBaseDepth(profile), 0),
    sigmaV: verticalEffectiveStress(profile, z), gammaSub: l.gammaSub, su, phiDeg: l.phiDeg, Dr: l.Dr,
    eps50: l.eps50, OCR: l.OCR, G0, Vs: l.Vs, k: l.k, PI: l.PI, St: l.St, layer: l,
  };
}

/**
 * Attach general and local scour to a profile (returns a new profile; depths stay referenced to the original
 * mudline). General scour removes overburden everywhere; local scour removes soil support around the pile and
 * reduces overburden linearly down to the overburden reduction depth.
 * Reference: API RP 2A-WSD 21st ed. C6.3 / API RP 2GEO (2011) scour guidance: for an isolated pile, in the absence
 * of project data, local scour depth 1.5 D and overburden reduction depth 6 D below the general scour level.
 * Validity: general + local < profile depth; reductionDepth >= local.
 * @param {object} profile profile from createProfile
 * @param {{ general?: number, local?: number, D?: number, reductionDepth?: number }} opts general [m] general
 *   scour depth (default 0); local [m] local scour depth (default 1.5 D); D [m] pile diameter (needed for the
 *   defaults); reductionDepth [m] overburden reduction depth below the general scour level (default 6 D)
 * @returns {object} new profile with scour = { general, local, reductionDepth } [m]
 */
export function applyScour(profile, opts = {}) {
  const general = opts.general ?? 0;
  if (!num(general) || general < 0) throw new RangeError('general: must be >= 0 m');
  if (opts.local === undefined && !num(opts.D)) throw new RangeError('D: required when local scour depth is not given');
  const local = opts.local ?? 1.5 * opts.D;
  if (!num(local) || local < 0) throw new RangeError('local: must be >= 0 m');
  if (opts.reductionDepth === undefined && !num(opts.D)) throw new RangeError('D: required when reductionDepth is not given');
  const reductionDepth = opts.reductionDepth ?? 6 * opts.D;
  if (!num(reductionDepth) || reductionDepth < local) throw new RangeError('reductionDepth: must be >= local scour depth');
  if (general + local >= profile.depth) throw new RangeError('general: scour removes the whole profile');
  return { ...profile, scour: { general, local, reductionDepth } };
}

/**
 * Sample σ′v and the main strength properties on a depth grid (for plotting).
 * Reference: see verticalEffectiveStress.
 * @param {object} profile profile
 * @param {{ dz?: number }} [opts] dz [m] sampling interval (default 0.5); layer boundaries are always included
 * @returns {{ z: number[], sigmaV: number[], su: number[], phiDeg: number[], type: string[] }} z [m], sigmaV [Pa],
 *   su [Pa] (NaN where undefined), phiDeg [deg] (NaN where undefined)
 */
export function sampleProfile(profile, opts = {}) {
  const dz = opts.dz ?? 0.5;
  if (!num(dz) || dz <= 0) throw new RangeError('dz: must be > 0 m');
  const out = { z: [], sigmaV: [], su: [], phiDeg: [], type: [] };
  const push = (z) => {
    const p = propertiesAt(profile, z);
    out.z.push(z); out.sigmaV.push(p.sigmaV); out.su.push(p.su ?? NaN); out.phiDeg.push(p.phiDeg ?? NaN); out.type.push(p.type);
  };
  for (const l of profile.layers) {
    const n = Math.max(1, Math.ceil(l.thickness / dz - 1e-9));
    // stop just short of the interface so the sample belongs to this layer, then the next layer starts on it
    for (let i = 0; i < n; i++) push(l.top + (l.thickness * i) / n);
    push(l.bottom - 1e-9 * l.thickness);
  }
  return out;
}

/**
 * Generic example profiles for demonstration only. They are NOT site data and are not taken from any project or
 * standard; values are round numbers inside commonly quoted ranges for each soil description.
 * Each entry: { label, layers } with layers as accepted by createProfile.
 */
export const EXAMPLE_PROFILES = Object.freeze({
  softClayNC: {
    label: 'Generic example — soft normally consolidated clay',
    layers: [{ type: 'clay', top: 0, thickness: 80, gammaSub: 6000, suTop: 2000, suGrad: 1500, eps50: 0.02, OCR: 1, PI: 40, St: 3 }],
  },
  stiffClayOC: {
    label: 'Generic example — stiff overconsolidated clay',
    layers: [{ type: 'clay', top: 0, thickness: 80, gammaSub: 9000, suTop: 100e3, suGrad: 2000, eps50: 0.005, OCR: 8, PI: 25, St: 2 }],
  },
  mediumDenseSand: {
    label: 'Generic example — medium-dense sand',
    layers: [{ type: 'sand', top: 0, thickness: 80, gammaSub: 9500, phiDeg: 33, Dr: 0.55, apiSand: 'medium-dense sand' }],
  },
  denseSand: {
    label: 'Generic example — dense sand',
    layers: [{ type: 'sand', top: 0, thickness: 80, gammaSub: 10000, phiDeg: 38, Dr: 0.8, apiSand: 'dense sand' }],
  },
  layeredSandClay: {
    label: 'Generic example — layered sand over clay over dense sand',
    layers: [
      { type: 'sand', top: 0, thickness: 8, gammaSub: 9500, phiDeg: 32, Dr: 0.5, apiSand: 'medium-dense sand' },
      { type: 'clay', top: 8, thickness: 17, gammaSub: 8000, suTop: 40e3, suGrad: 2000, eps50: 0.01, OCR: 2, PI: 30, St: 3 },
      { type: 'sand', top: 25, thickness: 55, gammaSub: 10000, phiDeg: 37, Dr: 0.8, apiSand: 'dense sand' },
    ],
  },
});

/**
 * Build one of the generic example profiles (demonstration only, see EXAMPLE_PROFILES).
 * @param {string} name key of EXAMPLE_PROFILES
 * @returns {object} profile from createProfile
 */
export function exampleProfile(name) {
  const e = EXAMPLE_PROFILES[name];
  if (!e) throw new RangeError(`name: unknown example profile '${name}'`);
  return createProfile(e.layers.map((l) => ({ ...l })));
}
