/**
 * Cross-section property calculators for frame analysis and steel code checks.
 *
 * Axis convention (shared with frame3d.js): x is the member axis; y and z are the section axes.
 * Iy is the second moment of area about the y axis (bending in the x–z plane, moment My) and is the
 * strong axis of I-sections whose web lies along z. Asy / Asz are the shear areas for shear forces
 * along y / z. All lengths in m, areas m², second moments m⁴, section moduli m³, warping constant m⁶.
 *
 * References: Roark's Formulas for Stress and Strain (8th ed.) Table A.1; Cowper (1966) "The shear
 * coefficient in Timoshenko's beam theory", J. Appl. Mech. 33; EN 1993-1-1 §6.2.6(3) for shear areas of
 * I-sections; Timoshenko & Gere, Theory of Elastic Stability, for warping constants.
 */

function pos(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) {
    throw new RangeError(`${name} must be a positive finite number`);
  }
  return v;
}

function nonNeg(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v >= 0)) {
    throw new RangeError(`${name} must be a non-negative finite number`);
  }
  return v;
}

function finish(p) {
  p.iy = Math.sqrt(p.Iy / p.A);
  p.iz = Math.sqrt(p.Iz / p.A);
  if (p.Ip === undefined) p.Ip = p.Iy + p.Iz;
  return p;
}

/**
 * Circular hollow section (tube).
 * Shear coefficient: Cowper (1966) for a hollow circle, κ = 6(1+ν)(1+m²)² / [(7+6ν)(1+m²)² + (20+12ν)m²], m = d/D.
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.t wall thickness [m], 0 < t <= D/2
 * @param {number} [o.nu=0.3] Poisson's ratio [-] (shear coefficient only)
 * @returns {{type:string,D:number,t:number,A:number,Iy:number,Iz:number,J:number,Ip:number,Asy:number,Asz:number,
 *   Wely:number,Welz:number,Wply:number,Wplz:number,iy:number,iz:number,Iw:number,cy:number,cz:number,perimeter:number}}
 *   A, Asy, Asz [m²]; Iy, Iz, J, Ip [m⁴]; Wel*, Wpl* [m³]; iy, iz, cy, cz (extreme fibre distances), perimeter [m]; Iw [m⁶]
 */
export function tubeSection({ D, t, nu = 0.3 }) {
  pos('D', D); pos('t', t);
  if (t > D / 2) throw new RangeError('t must not exceed D/2');
  const d = D - 2 * t;
  const A = Math.PI / 4 * (D * D - d * d);
  const I = Math.PI / 64 * (D ** 4 - d ** 4);
  const m = d / D;
  const kappa = 6 * (1 + nu) * (1 + m * m) ** 2 / ((7 + 6 * nu) * (1 + m * m) ** 2 + (20 + 12 * nu) * m * m);
  const Wpl = (D ** 3 - d ** 3) / 6;
  return finish({
    type: 'tube', D, t, A, Iy: I, Iz: I, J: 2 * I, Ip: 2 * I, Asy: kappa * A, Asz: kappa * A,
    Wely: 2 * I / D, Welz: 2 * I / D, Wply: Wpl, Wplz: Wpl, Iw: 0, cy: D / 2, cz: D / 2, perimeter: Math.PI * D,
  });
}

/**
 * Solid circular section. Shear coefficient κ = 6(1+ν)/(7+6ν) (Cowper 1966).
 * @param {object} o
 * @param {number} o.D diameter [m]
 * @param {number} [o.nu=0.3] Poisson's ratio [-]
 * @returns {object} same fields as tubeSection (SI)
 */
export function solidRoundSection({ D, nu = 0.3 }) {
  pos('D', D);
  const A = Math.PI * D * D / 4;
  const I = Math.PI * D ** 4 / 64;
  const kappa = 6 * (1 + nu) / (7 + 6 * nu);
  return finish({
    type: 'round', D, A, Iy: I, Iz: I, J: 2 * I, Ip: 2 * I, Asy: kappa * A, Asz: kappa * A,
    Wely: 2 * I / D, Welz: 2 * I / D, Wply: D ** 3 / 6, Wplz: D ** 3 / 6, Iw: 0, cy: D / 2, cz: D / 2,
    perimeter: Math.PI * D,
  });
}

/**
 * Welded box section with two webs (thickness tw, along z) and two flanges (thickness tf, along y).
 * Torsion constant by Bredt's formula on the wall mid-line, J = 4 Am² / ∮ds/t. Shear areas are the wall
 * areas parallel to the shear force (webs for z, flanges for y).
 * @param {object} o
 * @param {number} o.h overall depth along z [m]
 * @param {number} o.b overall width along y [m]
 * @param {number} o.tw web thickness [m]
 * @param {number} o.tf flange thickness [m]
 * @returns {object} section properties (SI) as tubeSection
 */
export function boxSection({ h, b, tw, tf }) {
  pos('h', h); pos('b', b); pos('tw', tw); pos('tf', tf);
  if (2 * tf >= h) throw new RangeError('tf must be smaller than h/2');
  if (2 * tw >= b) throw new RangeError('tw must be smaller than b/2');
  const hi = h - 2 * tf;
  const bi = b - 2 * tw;
  const A = b * h - bi * hi;
  const Iy = (b * h ** 3 - bi * hi ** 3) / 12;
  const Iz = (h * b ** 3 - hi * bi ** 3) / 12;
  const hm = h - tf;
  const bm = b - tw;
  const J = 4 * (hm * bm) ** 2 / (2 * bm / tf + 2 * hm / tw);
  return finish({
    type: 'box', h, b, tw, tf, A, Iy, Iz, J, Asz: 2 * tw * hi, Asy: 2 * tf * bi,
    Wely: 2 * Iy / h, Welz: 2 * Iz / b, Wply: (b * h * h - bi * hi * hi) / 4, Wplz: (h * b * b - hi * bi * bi) / 4,
    Iw: 0, cy: b / 2, cz: h / 2, perimeter: 2 * (b + h),
  });
}

/**
 * Rectangular hollow section of uniform wall thickness with sharp corners (boxSection with tw = tf = t).
 * @param {object} o
 * @param {number} o.h depth along z [m]
 * @param {number} o.b width along y [m]
 * @param {number} o.t wall thickness [m]
 * @returns {object} section properties (SI)
 */
export function rectHollowSection({ h, b, t }) {
  const s = boxSection({ h, b, tw: t, tf: t });
  s.type = 'rhs'; s.t = t;
  return s;
}

/**
 * Solid rectangle. Torsion constant from the series fit J = a b³ [1/3 − 0.21 (b/a)(1 − b⁴/(12 a⁴))], a >= b
 * (Roark Table 10.1 case 4). Shear coefficient 5/6.
 * @param {object} o
 * @param {number} o.h depth along z [m]
 * @param {number} o.b width along y [m]
 * @returns {object} section properties (SI)
 */
export function rectSection({ h, b }) {
  pos('h', h); pos('b', b);
  const a = Math.max(h, b);
  const c = Math.min(h, b);
  const J = a * c ** 3 * (1 / 3 - 0.21 * (c / a) * (1 - c ** 4 / (12 * a ** 4)));
  const A = b * h;
  return finish({
    type: 'rect', h, b, A, Iy: b * h ** 3 / 12, Iz: h * b ** 3 / 12, J, Asy: 5 / 6 * A, Asz: 5 / 6 * A,
    Wely: b * h * h / 6, Welz: h * b * b / 6, Wply: b * h * h / 4, Wplz: h * b * b / 4, Iw: 0, cy: b / 2, cz: h / 2,
    perimeter: 2 * (b + h),
  });
}

/* Spandrel fillet (square r×r minus quarter circle) : area, centroid distance from the two tangent faces and
   second moment about its own centroidal axis parallel to a face. */
function fillet(r) {
  if (r <= 0) return { A: 0, c: 0, I0: 0 };
  const A = r * r * (1 - Math.PI / 4);
  const c = r * (10 - 3 * Math.PI) / (3 * (4 - Math.PI));
  const aq = Math.PI * r * r / 4;
  const cq = 4 * r / (3 * Math.PI);
  // quarter circle about the face line, by parallel axis from its own centroid
  const Iq = Math.PI * r ** 4 / 16 - aq * cq * cq + aq * (r - cq) ** 2;
  const Iface = r ** 4 / 3 - Iq;
  return { A, c, I0: Iface - A * c * c };
}

/**
 * Doubly symmetric I/H section, rolled-like (root radius r) or welded (r = 0). Web along z, strong axis y.
 * Area, second moments and plastic moduli include the root fillets exactly. Torsion constant is the
 * thin-walled sum J = (2 b tf³ + (h − 2 tf) tw³)/3 (fillet contribution neglected, slightly conservative).
 * Warping constant Iw = Iz (h − tf)²/4. Shear area Asz per EN 1993-1-1 §6.2.6(3)(a) for rolled sections
 * (A − 2 b tf + (tw + 2r) tf, not less than hw tw) and hw tw for welded sections; Asy = (5/6)·2 b tf.
 * @param {object} o
 * @param {number} o.h overall depth [m]
 * @param {number} o.b flange width [m]
 * @param {number} o.tw web thickness [m]
 * @param {number} o.tf flange thickness [m]
 * @param {number} [o.r=0] root radius [m]
 * @param {boolean} [o.welded] fabrication route; default true when r = 0
 * @returns {object} { type, h, b, tw, tf, r, welded, A, Iy, Iz, J, Iw, Asy, Asz, Wely, Welz, Wply, Wplz, iy, iz, Ip, cy, cz, hw, perimeter } (SI)
 */
export function iSection({ h, b, tw, tf, r = 0, welded }) {
  pos('h', h); pos('b', b); pos('tw', tw); pos('tf', tf); nonNeg('r', r);
  if (2 * tf >= h) throw new RangeError('tf must be smaller than h/2');
  if (tw > b) throw new RangeError('tw must not exceed b');
  const hw = h - 2 * tf;
  if (2 * r > hw || 2 * r > b - tw) throw new RangeError('r is too large for the section');
  const f = fillet(r);
  const A = 2 * b * tf + hw * tw + 4 * f.A;
  const zf = hw / 2 - f.c;           // fillet centroid height above the centroid
  const yf = tw / 2 + f.c;           // fillet centroid offset from the web axis
  const Iy = (b * h ** 3 - (b - tw) * hw ** 3) / 12 + 4 * (f.I0 + f.A * zf * zf);
  const Iz = (2 * tf * b ** 3 + hw * tw ** 3) / 12 + 4 * (f.I0 + f.A * yf * yf);
  const J = (2 * b * tf ** 3 + hw * tw ** 3) / 3;
  const Iw = Iz * (h - tf) ** 2 / 4;
  const Wply = b * tf * (h - tf) + tw * hw * hw / 4 + 4 * f.A * zf;
  const Wplz = tf * b * b / 2 + hw * tw * tw / 4 + 4 * f.A * yf;
  const isWelded = welded === undefined ? r === 0 : !!welded;
  const Asz = isWelded ? hw * tw : Math.max(A - 2 * b * tf + (tw + 2 * r) * tf, hw * tw);
  return finish({
    type: 'I', h, b, tw, tf, r, welded: isWelded, A, Iy, Iz, J, Iw, Asy: 5 / 6 * 2 * b * tf, Asz,
    Wely: 2 * Iy / h, Welz: 2 * Iz / b, Wply, Wplz, cy: b / 2, cz: h / 2, hw,
    perimeter: 4 * b + 2 * h - 2 * tw - (8 - 2 * Math.PI) * r,
  });
}

/**
 * Parametric European / British style I-section helper. Generates an HE/IPE/UB/UC-like section from
 * (h, b, tw, tf); when the root radius is omitted it is taken as r = tw (order of magnitude of rolled
 * profiles — pass the tabulated radius when a specific profile is meant).
 * @param {object} o
 * @param {number} o.h depth [m]
 * @param {number} o.b width [m]
 * @param {number} o.tw web thickness [m]
 * @param {number} o.tf flange thickness [m]
 * @param {number} [o.r] root radius [m]
 * @returns {object} iSection result with welded = false
 */
export function rolledISection({ h, b, tw, tf, r }) {
  return iSection({ h, b, tw, tf, r: r === undefined ? tw : r, welded: false });
}

/**
 * Dimensions of a few common European rolled I/H profiles (EN 10365 / former Euronorm 19-57 and 53-62):
 * h, b, tw, tf, r in m. Use with rolledISection or iSectionByName.
 */
export const EURO_I_PROFILES = Object.freeze({
  'IPE 200': { h: 0.200, b: 0.100, tw: 0.0056, tf: 0.0085, r: 0.012 },
  'IPE 300': { h: 0.300, b: 0.150, tw: 0.0071, tf: 0.0107, r: 0.015 },
  'IPE 400': { h: 0.400, b: 0.180, tw: 0.0086, tf: 0.0135, r: 0.021 },
  'IPE 500': { h: 0.500, b: 0.200, tw: 0.0102, tf: 0.0160, r: 0.021 },
  'IPE 600': { h: 0.600, b: 0.220, tw: 0.0120, tf: 0.0190, r: 0.024 },
  'HEA 200': { h: 0.190, b: 0.200, tw: 0.0065, tf: 0.0100, r: 0.018 },
  'HEA 300': { h: 0.290, b: 0.300, tw: 0.0085, tf: 0.0140, r: 0.027 },
  'HEB 200': { h: 0.200, b: 0.200, tw: 0.0090, tf: 0.0150, r: 0.018 },
  'HEB 300': { h: 0.300, b: 0.300, tw: 0.0110, tf: 0.0190, r: 0.027 },
  'HEB 400': { h: 0.400, b: 0.300, tw: 0.0135, tf: 0.0240, r: 0.027 },
  'HEB 600': { h: 0.600, b: 0.300, tw: 0.0155, tf: 0.0300, r: 0.027 },
});

/**
 * Section properties of a profile in EURO_I_PROFILES.
 * @param {string} name e.g. 'IPE 300'
 * @returns {object} iSection result (SI) with `name`
 */
export function iSectionByName(name) {
  const p = EURO_I_PROFILES[name];
  if (!p) throw new RangeError(`name: unknown profile "${name}"`);
  const s = rolledISection(p);
  s.name = name;
  return s;
}

/* Properties of a set of axis-parallel rectangles {y, z (centre), by, hz}: centroid, Iy, Iz. */
function rectAssembly(parts) {
  let A = 0; let Sy = 0; let Sz = 0;
  for (const p of parts) { const a = p.by * p.hz; A += a; Sy += a * p.z; Sz += a * p.y; }
  const zc = Sy / A; const yc = Sz / A;
  let Iy = 0; let Iz = 0; let Iyz = 0;
  for (const p of parts) {
    const a = p.by * p.hz;
    Iy += p.by * p.hz ** 3 / 12 + a * (p.z - zc) ** 2;
    Iz += p.hz * p.by ** 3 / 12 + a * (p.y - yc) ** 2;
    Iyz += a * (p.y - yc) * (p.z - zc);
  }
  return { A, yc, zc, Iy, Iz, Iyz };
}

/* Plastic modulus about the horizontal equal-area axis of stacked rectangles (bending about y). */
function plasticModulusY(parts) {
  let A = 0; let zmin = Infinity; let zmax = -Infinity;
  for (const p of parts) { A += p.by * p.hz; zmin = Math.min(zmin, p.z - p.hz / 2); zmax = Math.max(zmax, p.z + p.hz / 2); }
  const areaBelow = (zz) => {
    let s = 0;
    for (const p of parts) {
      const lo = p.z - p.hz / 2;
      s += p.by * Math.min(Math.max(zz - lo, 0), p.hz);
    }
    return s;
  };
  let lo = zmin; let hi = zmax;
  for (let i = 0; i < 80; i++) { const mid = 0.5 * (lo + hi); if (areaBelow(mid) < A / 2) lo = mid; else hi = mid; }
  const zp = 0.5 * (lo + hi);
  let W = 0;
  for (const p of parts) {
    const a = p.z - p.hz / 2; const c = p.z + p.hz / 2;
    if (c <= zp) W += p.by * p.hz * (zp - p.z);
    else if (a >= zp) W += p.by * p.hz * (p.z - zp);
    else W += p.by * ((zp - a) ** 2 + (c - zp) ** 2) / 2;
  }
  return { Wpl: W, zp };
}

/**
 * Welded monosymmetric I-section / plate girder: top flange (bft × tft), web (hw × tw), bottom flange
 * (bfb × tfb). Web along z; z measured upward from the bottom fibre. Warping constant for unequal flanges
 * Iw = hs² I1 I2 / (I1 + I2) with I1, I2 the flange second moments about z and hs the distance between flange
 * centroids (Timoshenko & Gere §5.3); shear centre located at hs·I1/(I1+I2) above the bottom flange centroid
 * (I1 = top flange).
 * @param {object} o
 * @param {number} o.bft top flange width [m]
 * @param {number} o.tft top flange thickness [m]
 * @param {number} o.hw clear web depth [m]
 * @param {number} o.tw web thickness [m]
 * @param {number} [o.bfb=o.bft] bottom flange width [m]
 * @param {number} [o.tfb=o.tft] bottom flange thickness [m]
 * @returns {object} { type, h, A, zc (centroid above bottom fibre), zsc (shear centre above bottom fibre), Iy, Iz, J, Iw,
 *   Wely (minimum), WelyTop, WelyBot, Welz, Wply, Wplz, zp (plastic neutral axis above bottom fibre), Asy, Asz, iy, iz, webSlenderness (hw/tw) } (SI)
 */
export function weldedISection({ bft, tft, hw, tw, bfb = bft, tfb = tft }) {
  pos('bft', bft); pos('tft', tft); pos('hw', hw); pos('tw', tw); pos('bfb', bfb); pos('tfb', tfb);
  const h = tfb + hw + tft;
  const parts = [
    { y: 0, z: tfb / 2, by: bfb, hz: tfb },
    { y: 0, z: tfb + hw / 2, by: tw, hz: hw },
    { y: 0, z: tfb + hw + tft / 2, by: bft, hz: tft },
  ];
  const g = rectAssembly(parts);
  const pm = plasticModulusY(parts);
  const I1 = tft * bft ** 3 / 12;
  const I2 = tfb * bfb ** 3 / 12;
  const hs = h - tft / 2 - tfb / 2;
  const J = (bft * tft ** 3 + bfb * tfb ** 3 + hw * tw ** 3) / 3;
  const WelyTop = g.Iy / (h - g.zc);
  const WelyBot = g.Iy / g.zc;
  const bmax = Math.max(bft, bfb);
  return finish({
    type: 'weldedI', h, b: bmax, bft, tft, bfb, tfb, hw, tw, tf: Math.min(tft, tfb), welded: true, r: 0,
    A: g.A, zc: g.zc, zsc: tfb / 2 + hs * I1 / (I1 + I2),
    Iy: g.Iy, Iz: g.Iz, J, Iw: hs * hs * I1 * I2 / (I1 + I2), WelyTop, WelyBot, Wely: Math.min(WelyTop, WelyBot),
    Welz: 2 * g.Iz / bmax, Wply: pm.Wpl, zp: pm.zp, Wplz: (tft * bft * bft + tfb * bfb * bfb + hw * tw * tw) / 4,
    Asz: hw * tw, Asy: 5 / 6 * (bft * tft + bfb * tfb), cy: bmax / 2, cz: Math.max(g.zc, h - g.zc),
    webSlenderness: hw / tw, perimeter: 2 * (bft + bfb + h) - 2 * tw,
  });
}

/**
 * Plate girder with optional transverse stiffener spacing recorded for shear-buckling checks; geometry and
 * properties as weldedISection.
 * @param {object} o weldedISection arguments plus
 * @param {number} [o.stiffenerSpacing] transverse stiffener spacing a [m]
 * @returns {object} weldedISection result with type 'plateGirder' and panelAspect = a/hw when a is given
 */
export function plateGirderSection(o) {
  const s = weldedISection(o);
  s.type = 'plateGirder';
  if (o.stiffenerSpacing !== undefined) {
    pos('stiffenerSpacing', o.stiffenerSpacing);
    s.stiffenerSpacing = o.stiffenerSpacing;
    s.panelAspect = o.stiffenerSpacing / s.hw;
  }
  return s;
}

/**
 * Parallel-flange channel; web along z at the left, flanges pointing in +y. Thin-walled shear centre
 * e = 3 b'² tf / (6 b' tf + h' tw) measured from the web mid-line away from the flanges, and warping constant
 * Iw = tf b'³ h'² (3 b' tf + 2 h' tw) / [12 (6 b' tf + h' tw)], with b' = b − tw/2, h' = h − tf
 * (Timoshenko & Gere; Roark Table 10.2).
 * @param {object} o
 * @param {number} o.h depth [m]
 * @param {number} o.b flange width [m]
 * @param {number} o.tw web thickness [m]
 * @param {number} o.tf flange thickness [m]
 * @returns {object} { A, yc (centroid from web back face), Iy, Iz, J, Iw, eShear (shear centre from web mid-line), Wely, Welz (minimum),
 *   Wply, Asy, Asz, iy, iz } (SI)
 */
export function channelSection({ h, b, tw, tf }) {
  pos('h', h); pos('b', b); pos('tw', tw); pos('tf', tf);
  if (2 * tf >= h) throw new RangeError('tf must be smaller than h/2');
  if (tw >= b) throw new RangeError('tw must be smaller than b');
  const hw = h - 2 * tf;
  const parts = [
    { y: tw / 2, z: 0, by: tw, hz: hw },
    { y: b / 2, z: (h - tf) / 2, by: b, hz: tf },
    { y: b / 2, z: -(h - tf) / 2, by: b, hz: tf },
  ];
  const g = rectAssembly(parts);
  const bp = b - tw / 2;
  const hp = h - tf;
  const e = 3 * bp * bp * tf / (6 * bp * tf + hp * tw);
  const Iw = tf * bp ** 3 * hp * hp * (3 * bp * tf + 2 * hp * tw) / (12 * (6 * bp * tf + hp * tw));
  return finish({
    type: 'channel', h, b, tw, tf, A: g.A, yc: g.yc, Iy: g.Iy, Iz: g.Iz, J: (2 * b * tf ** 3 + hw * tw ** 3) / 3, Iw,
    eShear: e, Wely: 2 * g.Iy / h, Welz: g.Iz / Math.max(g.yc, b - g.yc), Wply: b * tf * (h - tf) + tw * hw * hw / 4,
    Asz: hw * tw, Asy: 5 / 6 * 2 * b * tf, cy: Math.max(g.yc, b - g.yc), cz: h / 2, perimeter: 4 * b + 2 * h - 2 * tw,
  });
}

/**
 * Angle (L) section with legs h (along z) and b (along y), uniform thickness t, sharp corners. Returns the
 * geometric-axis properties, product of inertia and principal second moments (Iu major, Iv minor).
 * @param {object} o
 * @param {number} o.h vertical leg [m]
 * @param {number} o.b horizontal leg [m]
 * @param {number} o.t thickness [m]
 * @returns {object} { A, yc, zc (centroid from heel), Iy, Iz, Iyz, Iu, Iv, alpha (principal angle from y to u) [rad], J, iy, iz, iu, iv,
 *   Wely, Welz (minimum elastic moduli about geometric axes) } (SI)
 */
export function angleSection({ h, b, t }) {
  pos('h', h); pos('b', b); pos('t', t);
  if (t >= Math.min(h, b)) throw new RangeError('t must be smaller than both legs');
  const parts = [
    { y: t / 2, z: h / 2, by: t, hz: h },
    { y: t + (b - t) / 2, z: t / 2, by: b - t, hz: t },
  ];
  const g = rectAssembly(parts);
  const avg = (g.Iy + g.Iz) / 2;
  const rad = Math.hypot((g.Iy - g.Iz) / 2, g.Iyz);
  const Iu = avg + rad;
  const Iv = avg - rad;
  const alpha = 0.5 * Math.atan2(-2 * g.Iyz, g.Iy - g.Iz);
  const s = finish({
    type: 'angle', h, b, t, A: g.A, yc: g.yc, zc: g.zc, Iy: g.Iy, Iz: g.Iz, Iyz: g.Iyz, Iu, Iv, alpha,
    J: (h + b - t) * t ** 3 / 3, Iw: 0, Wely: g.Iy / Math.max(g.zc, h - g.zc), Welz: g.Iz / Math.max(g.yc, b - g.yc),
    Asz: 5 / 6 * h * t, Asy: 5 / 6 * b * t, cy: Math.max(g.yc, b - g.yc), cz: Math.max(g.zc, h - g.zc),
    perimeter: 2 * (h + b),
  });
  s.iu = Math.sqrt(Iu / g.A);
  s.iv = Math.sqrt(Iv / g.A);
  return s;
}

/**
 * Tee section: flange b × tf on top, stem (h − tf) × tw below. z measured upward from the stem toe.
 * @param {object} o
 * @param {number} o.h overall depth [m]
 * @param {number} o.b flange width [m]
 * @param {number} o.tw stem thickness [m]
 * @param {number} o.tf flange thickness [m]
 * @returns {object} { A, zc (centroid above stem toe), Iy, Iz, J, WelyTop, WelyBot, Wely, Welz, Wply, zp, Wplz, Asy, Asz, iy, iz } (SI)
 */
export function teeSection({ h, b, tw, tf }) {
  pos('h', h); pos('b', b); pos('tw', tw); pos('tf', tf);
  if (tf >= h) throw new RangeError('tf must be smaller than h');
  const hs = h - tf;
  const parts = [
    { y: 0, z: hs / 2, by: tw, hz: hs },
    { y: 0, z: hs + tf / 2, by: b, hz: tf },
  ];
  const g = rectAssembly(parts);
  const pm = plasticModulusY(parts);
  const WelyTop = g.Iy / (h - g.zc);
  const WelyBot = g.Iy / g.zc;
  return finish({
    type: 'tee', h, b, tw, tf, A: g.A, zc: g.zc, Iy: g.Iy, Iz: g.Iz, J: (b * tf ** 3 + hs * tw ** 3) / 3, Iw: 0,
    WelyTop, WelyBot, Wely: Math.min(WelyTop, WelyBot), Welz: 2 * g.Iz / b, Wply: pm.Wpl, zp: pm.zp,
    Wplz: (tf * b * b + hs * tw * tw) / 4, Asz: 5 / 6 * hs * tw, Asy: 5 / 6 * b * tf, cy: b / 2,
    cz: Math.max(g.zc, h - g.zc), perimeter: 2 * (b + h),
  });
}

/**
 * Stiffened-plate strip: one stiffener (web hw × tw with optional flange bf × tf; bf = 0 gives a flat bar)
 * attached to a plate strip of width s (stiffener spacing or effective width) and thickness tp. Bending about
 * the axis parallel to the plate. z measured from the outer plate surface towards the stiffener.
 * @param {object} o
 * @param {number} o.s plate strip width [m]
 * @param {number} o.tp plate thickness [m]
 * @param {number} o.hw stiffener web depth [m]
 * @param {number} o.tw stiffener web thickness [m]
 * @param {number} [o.bf=0] stiffener flange width [m]
 * @param {number} [o.tf=0] stiffener flange thickness [m]
 * @returns {object} { A, Aplate, Astiff, zc (centroid from plate outer surface), Iy, Iz, WelPlate, WelStiff (elastic moduli at the plate
 *   surface and stiffener tip), Wply, zp, J, iy, iz, Asz } (SI)
 */
export function stiffenedPlateStrip({ s, tp, hw, tw, bf = 0, tf = 0 }) {
  pos('s', s); pos('tp', tp); pos('hw', hw); pos('tw', tw); nonNeg('bf', bf); nonNeg('tf', tf);
  const parts = [
    { y: 0, z: tp / 2, by: s, hz: tp },
    { y: 0, z: tp + hw / 2, by: tw, hz: hw },
  ];
  const hasFlange = bf > 0 && tf > 0;
  if (hasFlange) parts.push({ y: 0, z: tp + hw + tf / 2, by: bf, hz: tf });
  const h = tp + hw + (hasFlange ? tf : 0);
  const g = rectAssembly(parts);
  const pm = plasticModulusY(parts);
  return finish({
    type: 'stiffenedPlate', s, tp, hw, tw, bf, tf, h, A: g.A, Aplate: s * tp, Astiff: g.A - s * tp, zc: g.zc,
    Iy: g.Iy, Iz: g.Iz, WelPlate: g.Iy / g.zc, WelStiff: g.Iy / (h - g.zc), Wely: Math.min(g.Iy / g.zc, g.Iy / (h - g.zc)),
    Welz: 2 * g.Iz / s, Wply: pm.Wpl, zp: pm.zp, J: (s * tp ** 3 + hw * tw ** 3 + (hasFlange ? bf * tf ** 3 : 0)) / 3,
    Iw: 0, Asz: hw * tw, Asy: 5 / 6 * s * tp, cy: s / 2, cz: Math.max(g.zc, h - g.zc),
  });
}

/**
 * Area properties of a general polygonal section by Green's theorem (shoelace integrals). The outer loop
 * may be given in either orientation; holes are subtracted. No torsion constant or plastic modulus is
 * returned (they are not closed-form for a general polygon).
 * @param {object} o
 * @param {Array<[number,number]>|Float64Array} o.outer outer boundary vertices (y, z) [m]; flat array of 2n numbers accepted
 * @param {Array<Array<[number,number]>|Float64Array>} [o.holes] inner boundaries [m]
 * @returns {{A:number,yc:number,zc:number,Iy:number,Iz:number,Iyz:number,Iu:number,Iv:number,alpha:number,Ip:number,iy:number,iz:number,
 *   Wely:number,Welz:number,cy:number,cz:number,perimeter:number}} centroidal properties (SI); alpha = principal angle [rad]
 */
export function polygonSection({ outer, holes = [] }) {
  const toPts = (loop) => {
    if (!loop || loop.length === 0) throw new RangeError('outer must contain at least three vertices');
    if (typeof loop[0] === 'number') {
      const p = [];
      for (let i = 0; i + 1 < loop.length; i += 2) p.push([loop[i], loop[i + 1]]);
      return p;
    }
    return loop.map((q) => [q[0], q[1]]);
  };
  const integ = (pts) => {
    let A = 0; let Sy = 0; let Sz = 0; let Iyy = 0; let Izz = 0; let Iyz = 0; let per = 0;
    const n = pts.length;
    for (let i = 0; i < n; i++) {
      const [y0, z0] = pts[i];
      const [y1, z1] = pts[(i + 1) % n];
      const c = y0 * z1 - y1 * z0;
      A += c / 2;
      Sz += (y0 + y1) * c / 6;                      // ∫ y dA
      Sy += (z0 + z1) * c / 6;                      // ∫ z dA
      Iyy += (z0 * z0 + z0 * z1 + z1 * z1) * c / 12; // ∫ z² dA
      Izz += (y0 * y0 + y0 * y1 + y1 * y1) * c / 12; // ∫ y² dA
      Iyz += (y0 * z1 + 2 * y0 * z0 + 2 * y1 * z1 + y1 * z0) * c / 24;
      per += Math.hypot(y1 - y0, z1 - z0);
    }
    const s = A < 0 ? -1 : 1;
    return { A: s * A, Sy: s * Sy, Sz: s * Sz, Iyy: s * Iyy, Izz: s * Izz, Iyz: s * Iyz, per };
  };
  const po = toPts(outer);
  if (po.length < 3) throw new RangeError('outer must contain at least three vertices');
  const t = integ(po);
  for (const hl of holes) {
    const q = integ(toPts(hl));
    t.A -= q.A; t.Sy -= q.Sy; t.Sz -= q.Sz; t.Iyy -= q.Iyy; t.Izz -= q.Izz; t.Iyz -= q.Iyz;
  }
  if (!(t.A > 0)) throw new RangeError('outer encloses no net area');
  const yc = t.Sz / t.A;
  const zc = t.Sy / t.A;
  const Iy = t.Iyy - t.A * zc * zc;
  const Iz = t.Izz - t.A * yc * yc;
  const Iyz = t.Iyz - t.A * yc * zc;
  const avg = (Iy + Iz) / 2;
  const rad = Math.hypot((Iy - Iz) / 2, Iyz);
  let cy = 0; let cz = 0;
  for (const [y, z] of po) { cy = Math.max(cy, Math.abs(y - yc)); cz = Math.max(cz, Math.abs(z - zc)); }
  return finish({
    type: 'polygon', A: t.A, yc, zc, Iy, Iz, Iyz, Iu: avg + rad, Iv: avg - rad,
    alpha: 0.5 * Math.atan2(-2 * Iyz, Iy - Iz), Wely: Iy / cz, Welz: Iz / cy, cy, cz, perimeter: t.per,
  });
}

/**
 * Concrete- or grout-filled steel tube, optionally with an inner steel tube (grouted pile-in-leg):
 * transformed-section stiffness referred to the steel modulus. EA = Es As + kc Ec Ac, EI = Es Is + kc Ec Ic.
 * `concreteFactor` kc = 1 gives the uncracked transformed section; EN 1994-1-1 §6.7.3.3(3) uses 0.6 for the
 * effective flexural stiffness of composite columns.
 * @param {object} o
 * @param {number} o.D outer tube diameter [m]
 * @param {number} o.t outer tube thickness [m]
 * @param {number} o.Es steel modulus [Pa]
 * @param {number} o.Ec infill modulus [Pa]
 * @param {number} [o.concreteFactor=1] stiffness factor on the infill [-], 0..1
 * @param {number} [o.innerD] inner tube outer diameter [m]
 * @param {number} [o.innerT] inner tube thickness [m]
 * @param {boolean} [o.innerFilled=false] inner tube bore also filled
 * @param {number} [o.rhoS=7850] steel density [kg/m³]
 * @param {number} [o.rhoC=2400] infill density [kg/m³]
 * @returns {{EA:number,EI:number,A:number,Iy:number,Iz:number,J:number,Asy:number,Asz:number,As:number,Ac:number,Is:number,Ic:number,
 *   massPerLength:number,n:number}} EA [N], EI [N·m²]; A, Iy, Iz, J transformed to steel [m², m⁴]; steel and infill areas/inertias;
 *   mass per length [kg/m]; modular ratio n = Es/Ec
 */
export function compositeTubeSection({
  D, t, Es, Ec, concreteFactor = 1, innerD, innerT, innerFilled = false, rhoS = 7850, rhoC = 2400,
}) {
  pos('Es', Es); pos('Ec', Ec);
  if (!(concreteFactor >= 0 && concreteFactor <= 1)) throw new RangeError('concreteFactor must lie in [0, 1]');
  const outer = tubeSection({ D, t });
  const di = D - 2 * t;
  let As = outer.A; let Is = outer.Iy;
  let Ac = Math.PI * di * di / 4;
  let Ic = Math.PI * di ** 4 / 64;
  if (innerD !== undefined) {
    pos('innerD', innerD); pos('innerT', innerT);
    if (innerD >= di) throw new RangeError('innerD must be smaller than the bore of the outer tube');
    const inner = tubeSection({ D: innerD, t: innerT });
    As += inner.A; Is += inner.Iy;
    const bore = innerD - 2 * innerT;
    const Asolid = Math.PI * innerD * innerD / 4;
    const Isolid = Math.PI * innerD ** 4 / 64;
    Ac -= Asolid; Ic -= Isolid;
    if (innerFilled) { Ac += Math.PI * bore * bore / 4; Ic += Math.PI * bore ** 4 / 64; }
  }
  const n = Es / Ec;
  const A = As + concreteFactor * Ac / n;
  const I = Is + concreteFactor * Ic / n;
  return {
    type: 'compositeTube', D, t, EA: Es * A, EI: Es * I, A, Iy: I, Iz: I, J: 2 * Is, Ip: 2 * I, Asy: outer.Asy, Asz: outer.Asz,
    As, Ac, Is, Ic, n, massPerLength: rhoS * As + rhoC * Ac, iy: Math.sqrt(I / A), iz: Math.sqrt(I / A),
    Wely: 2 * I / D, Welz: 2 * I / D, cy: D / 2, cz: D / 2,
  };
}

/**
 * Compact library of tubular sizes commonly used for offshore jacket legs, braces, piles and conductors
 * (outer diameter × wall thickness, m). Sizes are nominal line-pipe / rolled-can dimensions.
 */
export const OFFSHORE_TUBULARS = Object.freeze([
  ['273x12.7', 0.2731, 0.0127], ['324x12.7', 0.3239, 0.0127], ['406x12.7', 0.4064, 0.0127], ['406x19.1', 0.4064, 0.0191],
  ['457x12.7', 0.4572, 0.0127], ['457x19.1', 0.4572, 0.0191], ['508x12.7', 0.5080, 0.0127], ['508x19.1', 0.5080, 0.0191],
  ['610x12.7', 0.6096, 0.0127], ['610x19.1', 0.6096, 0.0191], ['610x25.4', 0.6096, 0.0254], ['762x15.9', 0.7620, 0.0159],
  ['762x25.4', 0.7620, 0.0254], ['914x19.1', 0.9144, 0.0191], ['914x25.4', 0.9144, 0.0254], ['1067x25.4', 1.0668, 0.0254],
  ['1219x25.4', 1.2192, 0.0254], ['1219x38.1', 1.2192, 0.0381], ['1372x31.8', 1.3716, 0.0318], ['1524x38.1', 1.5240, 0.0381],
  ['1829x44.5', 1.8288, 0.0445], ['2134x50.8', 2.1336, 0.0508], ['2438x57.2', 2.4384, 0.0572],
].map(([name, D, t]) => Object.freeze({ name, D, t })));

/**
 * Section properties of a library tubular.
 * @param {string} name e.g. '914x25.4' (outer diameter × thickness in mm)
 * @returns {object} tubeSection result (SI) with `name`
 */
export function tubularByName(name) {
  const p = OFFSHORE_TUBULARS.find((q) => q.name === name);
  if (!p) throw new RangeError(`name: unknown tubular "${name}"`);
  const s = tubeSection({ D: p.D, t: p.t });
  s.name = name;
  return s;
}

/**
 * Extreme-fibre distances and elastic section moduli of a section property record, for stress recovery.
 * Uses, in order: the record's own moduli Wely / Welz; its fibre distances cz (for bending about y, measured
 * along z) and cy (for bending about z); otherwise — a generic record {A, Iy, Iz, J} as produced by a frame
 * import, which carries no shape — the distances are taken from the radii of gyration as c = √2·i,
 * cz = √(2 Iy/A), cy = √(2 Iz/A). That is exact for a thin-walled circular tube (i = R/√2) and is flagged by
 * `assumed`: it overstates the stress of an I-section about its strong axis (c ≈ 1.2 i) by about 15 % and
 * understates that of a solid rectangle (c = √3·i) by 18 %. Supply cy and cz when the shape is known.
 * Reference: elastic beam theory, σ = M c / I (e.g. Timoshenko & Gere, Mechanics of Materials).
 * @param {object} sec section property record with A [m²], Iy, Iz [m⁴] and optionally Wely, Welz [m³], cy, cz [m]
 * @returns {{cy:number, cz:number, Wely:number, Welz:number, assumed:boolean}} fibre distances [m], elastic
 *   moduli [m³]; assumed = at least one distance came from the radius of gyration
 */
export function fibreDistances(sec) {
  if (!sec || !(sec.A > 0) || !(sec.Iy > 0) || !(sec.Iz > 0)) throw new RangeError('section must have positive A, Iy and Iz');
  for (const k of ['cy', 'cz']) {
    if (sec[k] !== undefined && !(sec[k] > 0 && Number.isFinite(sec[k]))) throw new RangeError(`section.${k} must be a positive finite number`);
  }
  let assumed = false;
  let Wely; let cz;
  if (sec.Wely > 0) { Wely = sec.Wely; cz = sec.cz > 0 ? sec.cz : sec.Iy / Wely; }
  else {
    if (sec.cz > 0) cz = sec.cz; else { cz = Math.sqrt(2 * sec.Iy / sec.A); assumed = true; }
    Wely = sec.Iy / cz;
  }
  let Welz; let cy;
  if (sec.Welz > 0) { Welz = sec.Welz; cy = sec.cy > 0 ? sec.cy : sec.Iz / Welz; }
  else {
    if (sec.cy > 0) cy = sec.cy; else { cy = Math.sqrt(2 * sec.Iz / sec.A); assumed = true; }
    Welz = sec.Iz / cy;
  }
  return { cy, cz, Wely, Welz, assumed };
}

/**
 * Resolve a section description into a property record usable by frame3d: accepts a property record
 * ({A, Iy, Iz, J, ...}) unchanged, or a typed description ({type:'tube', D, t}, {type:'I', h, b, tw, tf, r},
 * {type:'box'|'rhs'|'rect'|'round'|'channel'|'angle'|'tee'|'weldedI', ...}).
 * @param {object} s section description
 * @returns {object} section property record (SI)
 */
export function resolveSection(s) {
  if (!s || typeof s !== 'object') throw new RangeError('section must be an object');
  if (typeof s.A === 'number' && typeof s.Iy === 'number' && typeof s.Iz === 'number') {
    if (s.J === undefined) throw new RangeError('section.J is required when properties are given directly');
    return s;
  }
  switch (s.type) {
    case 'tube': return tubeSection(s);
    case 'round': return solidRoundSection(s);
    case 'I': return iSection(s);
    case 'box': return boxSection(s);
    case 'rhs': return rectHollowSection(s);
    case 'rect': return rectSection(s);
    case 'channel': return channelSection(s);
    case 'angle': return angleSection(s);
    case 'tee': return teeSection(s);
    case 'weldedI': return weldedISection(s);
    default: throw new RangeError(`section.type "${s.type}" is not recognised`);
  }
}
