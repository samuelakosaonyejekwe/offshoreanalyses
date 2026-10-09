// Cross-section properties from dimensions, and a resolver from section designations (as they appear in SDNF,
// CIS/2, IFC and solver decks) to the contract's section records.
//
// Axis convention (Eurocode 3): x along the member, y–y the major (strong) principal axis, z–z the minor.
// So for an I-section Iy is bending about the axis parallel to the flanges. All values SI (m, m², m⁴).
//
// The catalogue is deliberately small: IPE 80–600, HE A / HE B 100–1000 (EN 10365 nominal dimensions) and a
// selection of AISC W shapes (ASTM A6 nominal dimensions, US customary). HE M, UB/UC, metric W designations,
// channels, angles and tees are NOT tabulated — such names are returned unresolved (name kept, no properties).

const MM = 1e-3, INCH = 0.0254;

/** Circular hollow section. @param {number} D outer diameter, m @param {number} t wall thickness, m */
export function tubeProps(D, t) {
  const d = Math.max(D - 2 * t, 0);
  const I = (Math.PI / 64) * (D ** 4 - d ** 4);
  return { A: (Math.PI / 4) * (D * D - d * d), Iy: I, Iz: I, J: 2 * I };
}

/** Solid circle. @param {number} D diameter, m */
export function circleProps(D) {
  const I = (Math.PI / 64) * D ** 4;
  return { A: (Math.PI / 4) * D * D, Iy: I, Iz: I, J: 2 * I };
}

/**
 * Solid rectangle, width b (along y) × depth h (along z); Iy is about the y axis (b·h³/12).
 * Torsion constant from the series approximation J = a·c³·(1/3 − 0.21·(c/a)·(1 − c⁴/(12a⁴))), a ≥ c
 * (Roark's Formulas for Stress and Strain, Table 10.1 case 4).
 */
export function rectProps(b, h) {
  const a = Math.max(b, h), c = Math.min(b, h);
  const J = a * c ** 3 * (1 / 3 - 0.21 * (c / a) * (1 - c ** 4 / (12 * a ** 4)));
  return { A: b * h, Iy: (b * h ** 3) / 12, Iz: (h * b ** 3) / 12, J };
}

/**
 * Rectangular hollow section, outer width b × depth h, wall tw (webs, the sides of depth h) and tf (flanges).
 * Sharp corners (corner radii neglected). Torsion constant by Bredt's formula on the wall mid-line.
 */
export function boxProps(b, h, tw, tf = tw) {
  const bi = Math.max(b - 2 * tw, 0), hi = Math.max(h - 2 * tf, 0);
  const Am = (b - tw) * (h - tf);
  return {
    A: b * h - bi * hi,
    Iy: (b * h ** 3 - bi * hi ** 3) / 12,
    Iz: (h * b ** 3 - hi * bi ** 3) / 12,
    J: (4 * Am * Am) / ((2 * (b - tw)) / tf + (2 * (h - tf)) / tw),
  };
}

/**
 * Doubly symmetric I / H section: depth h, flange width b, web tw, flange tf, root radius r (0 to neglect).
 * Area and second moments include the four root fillets; J is the thin-walled open-section sum
 * (2·b·tf³ + (h − tf)·tw³)/3, which neglects the fillet contribution (low by 10–30 % for rolled shapes).
 */
export function iSectionProps({ h, b, tw, tf, r = 0 }) {
  const hw = h - 2 * tf;
  const fa = (1 - Math.PI / 4) * r * r;               // one fillet area
  const fc = r * (1 - 2 / (3 * (4 - Math.PI)));        // fillet centroid distance from the web face / flange face
  const zc = hw / 2 - fc, yc = tw / 2 + fc;
  const fi = (1 / 3 - Math.PI / 16 - 1 / (9 * (4 - Math.PI))) * r ** 4; // fillet centroidal second moment
  return {
    A: 2 * b * tf + hw * tw + 4 * fa,
    Iy: (b * h ** 3 - (b - tw) * hw ** 3) / 12 + 4 * (fi + fa * zc * zc),
    Iz: (2 * tf * b ** 3 + hw * tw ** 3) / 12 + 4 * (fi + fa * yc * yc),
    J: (2 * b * tf ** 3 + (h - tf) * tw ** 3) / 3,
  };
}

// EN 10365 nominal dimensions, mm: [h, b, tw, tf, r]
const IPE = {
  80: [80, 46, 3.8, 5.2, 5], 100: [100, 55, 4.1, 5.7, 7], 120: [120, 64, 4.4, 6.3, 7], 140: [140, 73, 4.7, 6.9, 7],
  160: [160, 82, 5.0, 7.4, 9], 180: [180, 91, 5.3, 8.0, 9], 200: [200, 100, 5.6, 8.5, 12], 220: [220, 110, 5.9, 9.2, 12],
  240: [240, 120, 6.2, 9.8, 15], 270: [270, 135, 6.6, 10.2, 15], 300: [300, 150, 7.1, 10.7, 15], 330: [330, 160, 7.5, 11.5, 18],
  360: [360, 170, 8.0, 12.7, 18], 400: [400, 180, 8.6, 13.5, 21], 450: [450, 190, 9.4, 14.6, 21], 500: [500, 200, 10.2, 16.0, 21],
  550: [550, 210, 11.1, 17.2, 24], 600: [600, 220, 12.0, 19.0, 24],
};
const HEA = {
  100: [96, 100, 5, 8, 12], 120: [114, 120, 5, 8, 12], 140: [133, 140, 5.5, 8.5, 12], 160: [152, 160, 6, 9, 15],
  180: [171, 180, 6, 9.5, 15], 200: [190, 200, 6.5, 10, 18], 220: [210, 220, 7, 11, 18], 240: [230, 240, 7.5, 12, 21],
  260: [250, 260, 7.5, 12.5, 24], 280: [270, 280, 8, 13, 24], 300: [290, 300, 8.5, 14, 27], 320: [310, 300, 9, 15.5, 27],
  340: [330, 300, 9.5, 16.5, 27], 360: [350, 300, 10, 17.5, 27], 400: [390, 300, 11, 19, 27], 450: [440, 300, 11.5, 21, 27],
  500: [490, 300, 12, 23, 27], 550: [540, 300, 12.5, 24, 27], 600: [590, 300, 13, 25, 27], 650: [640, 300, 13.5, 26, 27],
  700: [690, 300, 14.5, 27, 27], 800: [790, 300, 15, 28, 30], 900: [890, 300, 16, 30, 30], 1000: [990, 300, 16.5, 31, 30],
};
const HEB = {
  100: [100, 100, 6, 10, 12], 120: [120, 120, 6.5, 11, 12], 140: [140, 140, 7, 12, 12], 160: [160, 160, 8, 13, 15],
  180: [180, 180, 8.5, 14, 15], 200: [200, 200, 9, 15, 18], 220: [220, 220, 9.5, 16, 18], 240: [240, 240, 10, 17, 21],
  260: [260, 260, 10, 17.5, 24], 280: [280, 280, 10.5, 18, 24], 300: [300, 300, 11, 19, 27], 320: [320, 300, 11.5, 20.5, 27],
  340: [340, 300, 12, 21.5, 27], 360: [360, 300, 12.5, 22.5, 27], 400: [400, 300, 13.5, 24, 27], 450: [450, 300, 14, 26, 27],
  500: [500, 300, 14.5, 28, 27], 550: [550, 300, 15, 29, 27], 600: [600, 300, 15.5, 30, 27], 650: [650, 300, 16, 31, 27],
  700: [700, 300, 17, 32, 27], 800: [800, 300, 17.5, 33, 30], 900: [900, 300, 18.5, 35, 30], 1000: [1000, 300, 19, 36, 30],
};

// AISC W shapes (ASTM A6), inches: 'depth x weight' → [d, bf, tw, tf]. The area is taken from the nominal
// weight (lb/ft ÷ 3.403, i.e. 490 lb/ft³ steel) because the table carries no fillet dimension; second moments
// are computed from the plates with the fillets neglected (1–3 % low).
export const AISC_W = {
  '6x15': [5.99, 5.99, 0.230, 0.260], '6x25': [6.38, 6.08, 0.320, 0.455],
  '8x18': [8.14, 5.25, 0.230, 0.330], '8x24': [7.93, 6.50, 0.245, 0.400], '8x31': [8.00, 8.00, 0.285, 0.435],
  '10x22': [10.2, 5.75, 0.240, 0.360], '10x33': [9.73, 7.96, 0.290, 0.435], '10x45': [10.1, 8.02, 0.350, 0.620], '10x49': [10.0, 10.0, 0.340, 0.560],
  '12x19': [12.2, 4.01, 0.235, 0.350], '12x26': [12.2, 6.49, 0.230, 0.380], '12x40': [11.9, 8.01, 0.295, 0.515], '12x53': [12.1, 10.0, 0.345, 0.575],
  '12x65': [12.1, 12.0, 0.390, 0.605], '12x72': [12.3, 12.0, 0.430, 0.670], '12x96': [12.7, 12.2, 0.550, 0.900],
  '14x22': [13.7, 5.00, 0.230, 0.335], '14x48': [13.8, 8.03, 0.340, 0.595], '14x61': [13.9, 10.0, 0.375, 0.645], '14x82': [14.3, 10.1, 0.510, 0.855],
  '14x90': [14.0, 14.5, 0.440, 0.710], '14x120': [14.5, 14.7, 0.590, 0.940], '14x132': [14.7, 14.7, 0.645, 1.03], '14x145': [14.8, 15.5, 0.680, 1.09],
  '14x176': [15.2, 15.7, 0.830, 1.31], '14x211': [15.7, 15.8, 0.980, 1.56], '14x257': [16.4, 16.0, 1.18, 1.89],
  '16x26': [15.7, 5.50, 0.250, 0.345], '16x31': [15.9, 5.53, 0.275, 0.440], '16x36': [15.9, 6.99, 0.295, 0.430], '16x40': [16.0, 7.00, 0.305, 0.505],
  '18x35': [17.7, 6.00, 0.300, 0.425], '18x40': [17.9, 6.02, 0.315, 0.525], '18x50': [18.0, 7.50, 0.355, 0.570], '18x60': [18.2, 7.56, 0.415, 0.695],
  '21x44': [20.7, 6.50, 0.350, 0.450], '21x50': [20.8, 6.53, 0.380, 0.535], '21x57': [21.1, 6.56, 0.405, 0.650], '21x62': [21.0, 8.24, 0.400, 0.615],
  '24x55': [23.6, 7.01, 0.395, 0.505], '24x62': [23.7, 7.04, 0.430, 0.590], '24x68': [23.7, 8.97, 0.415, 0.585], '24x76': [23.9, 8.99, 0.440, 0.680],
  '24x84': [24.1, 9.02, 0.470, 0.770], '24x104': [24.1, 12.8, 0.500, 0.750],
  '27x84': [26.7, 10.0, 0.460, 0.640], '27x94': [26.9, 10.0, 0.490, 0.745],
  '30x99': [29.7, 10.5, 0.520, 0.670], '30x108': [29.8, 10.5, 0.545, 0.760], '30x116': [30.0, 10.5, 0.565, 0.850],
  '33x118': [32.9, 11.5, 0.550, 0.740], '33x130': [33.1, 11.5, 0.580, 0.855],
  '36x135': [35.6, 12.0, 0.600, 0.790], '36x150': [35.9, 12.0, 0.625, 0.940], '36x160': [36.0, 12.0, 0.650, 1.02],
};

export const EURO_I = { IPE, HEA, HEB };

function iRecord(dims, scale, standard, designation) {
  const [h, b, tw, tf, r = 0] = dims.map((v) => v * scale);
  const p = iSectionProps({ h, b, tw, tf, r });
  return { ...p, shape: 'I', h, b, tw, tf, r, designation, standard };
}

const num = '(\\d+(?:[.,]\\d+)?)';
const SEP = '\\s*[xX×*/]\\s*';

/**
 * Resolve a section designation to properties.
 *
 * Recognised (case-insensitive, spaces optional):
 *  - circular hollow: `PIPE 610x19`, `CHS 508x16`, `TUBE 610x19`, `TUB610X19`, `PIP 610x19`, `Ø610x19`,
 *    `O 610x19`, `D610T19`, `PIPE 24"x0.5"` → `{ type: 'tube', D, t }`
 *  - rectangular / square hollow: `RHS 200x100x8`, `SHS 150x150x6` or `SHS 150x6`, `BOX 400x300x12`,
 *    `BOX 400x300x12x16` (h × b × tw × tf)
 *  - European rolled I/H: `IPE 300`, `HEA 200`, `HE 200 A`, `HEB300`, `HE300B`
 *  - AISC wide flange: `W14X90`, `W 14 x 90`
 *  - plates and bars: `PL 20x300`, `FLAT 300x20`, `FB 100x10`, `RECT 300x20`; `RD 50`, `ROD 50`, `BAR 50`
 *
 * Bare numbers are millimetres unless `opts.dimUnit` says otherwise ('mm' | 'cm' | 'm' | 'in'); a trailing
 * inch mark (") or `in` on a number forces inches for that designation.
 *
 * @param {string} name
 * @param {object} [opts] `{ dimUnit?: 'mm'|'cm'|'m'|'in' }`
 * @returns {{ matched: boolean, name: string, section: object|null, shape?: string, standard?: string }}
 *   `section` is `{ type: 'tube', D, t }` for circular hollow sections, otherwise
 *   `{ A, Iy, Iz, J, shape, ...dimensions }`; null when the designation is not recognised.
 */
export function resolveSection(name, opts = {}) {
  const raw = String(name ?? '').trim();
  const miss = { matched: false, name: raw, section: null };
  if (!raw) return miss;
  const inches = /\d\s*("|''|in\b|inch)/i.test(raw);
  const unit = inches ? INCH : ({ mm: MM, cm: 1e-2, m: 1, in: INCH })[opts.dimUnit || 'mm'] ?? MM;
  const s = raw.toUpperCase().replace(/"|''|\bINCH(ES)?\b|\bIN\b|\bMM\b/g, '').replace(/,(?=\d)/g, '.').replace(/\s+/g, ' ').trim();
  const f = (v) => parseFloat(v) * unit;
  let m;

  if ((m = new RegExp(`^(?:PIPE|PIP|CHS|TUBE|TUB|TUBULAR|CHHF|CHCF|HSS ?ROUND|Ø|⌀|O|P|PIPES?)[ _-]*${num}${SEP}${num}$`).exec(s)) ||
      (m = new RegExp(`^(?:OD|D)[ _-]*${num}[ _-]*(?:WT|T)[ _-]*${num}$`).exec(s))) {
    const D = f(m[1]), t = f(m[2]);
    if (t > 0 && 2 * t <= D) return { matched: true, name: raw, shape: 'tube', standard: 'dimensions', section: { type: 'tube', D, t } };
    return miss;
  }
  if ((m = new RegExp(`^(?:RHS|SHS|BOX|RHHF|RHCF|SHHF|SHCF|HSS)[ _-]*${num}${SEP}${num}(?:${SEP}${num})?(?:${SEP}${num})?$`).exec(s))) {
    let h, b, tw, tf;
    if (m[4] !== undefined) { h = f(m[1]); b = f(m[2]); tw = f(m[3]); tf = f(m[4]); } else if (m[3] !== undefined) { h = f(m[1]); b = f(m[2]); tw = tf = f(m[3]); } else { h = b = f(m[1]); tw = tf = f(m[2]); }
    if (!(tw > 0 && tf > 0 && 2 * tw < b && 2 * tf < h)) return miss;
    return { matched: true, name: raw, shape: 'box', standard: 'dimensions', section: { ...boxProps(b, h, tw, tf), shape: 'box', h, b, tw, tf } };
  }
  if ((m = /^(IPE|HEA|HEB)[ _-]*(\d{2,4})$/.exec(s)) || (m = /^HE[ _-]*(\d{2,4})[ _-]*([AB])$/.exec(s))) {
    const series = m[1] === 'IPE' || m[1] === 'HEA' || m[1] === 'HEB' ? m[1] : 'HE' + m[2];
    const size = m[1] === series ? m[2] : m[1];
    const dims = EURO_I[series] && EURO_I[series][Number(size)];
    if (!dims) return miss;
    const designation = series === 'IPE' ? `IPE ${Number(size)}` : `HE ${Number(size)} ${series[2]}`;
    return { matched: true, name: raw, shape: 'I', standard: 'EN 10365', section: iRecord(dims, MM, 'EN 10365', designation) };
  }
  if ((m = /^W[ _-]*(\d{1,2})[ _-]*[X×*][ _-]*(\d{1,3})$/.exec(s))) {
    const key = `${Number(m[1])}x${Number(m[2])}`;
    const dims = AISC_W[key];
    if (!dims) return miss;
    const rec = iRecord(dims, INCH, 'ASTM A6 / AISC', `W${key.toUpperCase().replace('X', 'X')}`);
    rec.A = (Number(m[2]) / 3.403) * INCH * INCH;
    rec.designation = `W${Number(m[1])}X${Number(m[2])}`;
    return { matched: true, name: raw, shape: 'I', standard: 'ASTM A6 / AISC', section: rec };
  }
  if ((m = new RegExp(`^(?:PL|PLT|PLATE|FLAT|FB|FL|RECT|BAR)[ _-]*${num}${SEP}${num}$`).exec(s))) {
    const a = f(m[1]), c = f(m[2]);
    const b = Math.min(a, c), h = Math.max(a, c);
    return { matched: true, name: raw, shape: 'rect', standard: 'dimensions', section: { ...rectProps(b, h), shape: 'rect', b, h } };
  }
  if ((m = new RegExp(`^(?:RD|ROD|RND|ROUND|BAR|SOLID ?ROUND)[ _-]*${num}$`).exec(s))) {
    const D = f(m[1]);
    return { matched: true, name: raw, shape: 'circle', standard: 'dimensions', section: { ...circleProps(D), shape: 'circle', D } };
  }
  return miss;
}
