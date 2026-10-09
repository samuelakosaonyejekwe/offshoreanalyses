// Time-dependent degradation: corrosion wastage, cathodic protection design, marine growth, wear and generic
// deterioration models. SI throughout; times are in years where the argument or field name says so.

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function positive(v, name) {
  need(Number.isFinite(v) && v > 0, `${name} must be a positive finite number`);
}

// ---------------------------------------------------------------------------------------------------------------
// Corrosion wastage
// ---------------------------------------------------------------------------------------------------------------

/**
 * Corrosion wastage depth at time t for one exposed surface.
 *  'linear'  d = r·(t − Tc) after the coating life Tc (constant rate r);
 *  'guedes-soares-garbatov'  d = d∞·[1 − exp(−(t − τc)/τt)] — Guedes Soares & Garbatov (1999) Marine Structures
 *      12, 425–445: long-term wastage d∞, coating life τc, transition time τt;
 *  'paik'  d = C1·(t − Tc)^C2 — Paik, Kim & Lee (1998) / Paik & Thayamballi (2003): C2 typically 0.3–1.5, C1 the
 *      wastage after one year of exposure.
 * Rates depend on zone and temperature; typical uncoated carbon steel values are 0.1–0.2 mm/yr immersed and
 * 0.3–0.5 mm/yr in the splash zone (order-of-magnitude literature values, not code values) — supply project values.
 * @param {object} o
 * @param {number} o.time age [yr]
 * @param {string} [o.model='linear']
 * @param {number} [o.coatingLife=0] Tc or τc [yr]
 * @param {number} [o.rate] constant rate r [m/yr] (linear)
 * @param {number} [o.dInf] long-term wastage d∞ [m] (Guedes Soares–Garbatov)
 * @param {number} [o.transitionTime] τt [yr] (Guedes Soares–Garbatov)
 * @param {number} [o.c1] C1 [m/yr^C2] (Paik)
 * @param {number} [o.c2=1] C2 [-] (Paik)
 * @returns {{depth:number, rate:number}} wastage [m] and instantaneous rate [m/yr]
 */
export function corrosionWastage({ time, model = 'linear', coatingLife = 0, rate, dInf, transitionTime, c1, c2 = 1 }) {
  need(time >= 0, 'time must be non-negative');
  need(coatingLife >= 0, 'coatingLife must be non-negative');
  const te = time - coatingLife;
  if (model === 'linear') {
    need(rate >= 0, 'rate must be non-negative');
    return te <= 0 ? { depth: 0, rate: 0 } : { depth: rate * te, rate };
  }
  if (model === 'guedes-soares-garbatov') {
    positive(dInf, 'dInf');
    positive(transitionTime, 'transitionTime');
    if (te <= 0) return { depth: 0, rate: 0 };
    const e = Math.exp(-te / transitionTime);
    return { depth: dInf * (1 - e), rate: (dInf / transitionTime) * e };
  }
  if (model === 'paik') {
    positive(c1, 'c1');
    positive(c2, 'c2');
    if (te <= 0) return { depth: 0, rate: 0 };
    return { depth: c1 * te ** c2, rate: c1 * c2 * te ** (c2 - 1) };
  }
  throw new RangeError(`model: unknown corrosion model "${model}"`);
}

/**
 * Remaining wall thickness and remaining life against a minimum-thickness criterion.
 * The wastage model applies to each corroding surface; `surfaces` = 2 for plating exposed on both sides.
 * @param {object} o
 * @param {number} o.thickness as-built thickness [m]
 * @param {number} o.minThickness minimum acceptable thickness [m]
 * @param {number} o.age current age [yr]
 * @param {number} [o.surfaces=1] number of corroding surfaces [-]
 * @param {object} o.corrosion arguments of {@link corrosionWastage} without `time`
 * @param {number} [o.horizon=200] search limit for the end of life [yr]
 * @returns {{remainingThickness:number, wastage:number, endOfLifeAge:number, remainingLife:number, ok:boolean}}
 *   thickness [m]; ages [yr]; endOfLifeAge = Infinity when the criterion is never reached within the horizon
 */
export function remainingWall({ thickness, minThickness, age, surfaces = 1, corrosion, horizon = 200 }) {
  positive(thickness, 'thickness');
  need(minThickness >= 0 && minThickness < thickness, 'minThickness must be in [0, thickness)');
  need(age >= 0, 'age must be non-negative');
  need(surfaces > 0, 'surfaces must be positive');
  const loss = (t) => surfaces * corrosionWastage({ ...corrosion, time: t }).depth;
  const allow = thickness - minThickness;
  const wastage = loss(age);
  let end = Infinity;
  if (loss(horizon) >= allow) {
    let lo = 0;
    let hi = horizon;
    for (let i = 0; i < 100; i++) {
      const mid = 0.5 * (lo + hi);
      if (loss(mid) >= allow) hi = mid; else lo = mid;
    }
    end = hi;
  }
  return { remainingThickness: thickness - wastage, wastage, endOfLifeAge: end, remainingLife: Math.max(0, end - age), ok: wastage < allow };
}

// ---------------------------------------------------------------------------------------------------------------
// Cathodic protection (galvanic anodes)
// ---------------------------------------------------------------------------------------------------------------

/**
 * Coating breakdown factors fc = a + b·t: initial (t = 0), mean (t = tf/2) and final (t = tf), capped at 1.
 * Reference: DNV-RP-B401 (coating breakdown factors for CP design). Constants per coating category (a, b for depth 0–30 m / > 30 m):
 * I: 0.10, 0.10 / 0.05;  II: 0.05, 0.025 / 0.015;  III: 0.02, 0.012 / 0.008. Bare steel: a = 1, b = 0.
 * Use {@link coatingBreakdownConstants} for these presets or pass project values.
 * @param {object} o
 * @param {number} o.a constant a [-]
 * @param {number} o.b constant b [1/yr]
 * @param {number} o.designLife tf [yr]
 * @returns {{initial:number, mean:number, final:number}} breakdown factors [-]
 */
export function coatingBreakdown({ a, b, designLife }) {
  need(a >= 0 && a <= 1, 'a must be in [0, 1]');
  need(b >= 0, 'b must be non-negative');
  need(designLife >= 0, 'designLife must be non-negative');
  const cap = (x) => Math.min(1, x);
  const final = a + b * designLife;
  // When the factor reaches 1 before the end of life the mean follows the capped line (DNV-RP-B401).
  let mean = a + (b * designLife) / 2;
  if (final > 1 && b > 0) mean = 1 - (1 - a) ** 2 / (2 * b * designLife);
  return { initial: cap(a), mean: cap(mean), final: cap(final) };
}

/**
 * Coating breakdown constants of DNV-RP-B401 by coating category and depth zone.
 * @param {string} category 'I' | 'II' | 'III' | 'bare'
 * @param {number} [depth=0] water depth of the surface [m]
 * @returns {{a:number, b:number}}
 */
export function coatingBreakdownConstants(category, depth = 0) {
  const deep = depth > 30;
  switch (category) {
    case 'I': return { a: 0.10, b: deep ? 0.05 : 0.10 };
    case 'II': return { a: 0.05, b: deep ? 0.015 : 0.025 };
    case 'III': return { a: 0.02, b: deep ? 0.008 : 0.012 };
    case 'bare': return { a: 1, b: 0 };
    default: throw new RangeError(`category: unknown coating category "${category}"`);
  }
}

/**
 * Anode resistance to remote seawater.
 *  'long-slender-standoff' (L >= 4r):  Ra = ρ/(2πL)·[ln(4L/r) − 1]
 *  'short-slender-standoff' (L < 4r):  Ra = ρ/(2πL)·[ln{(2L/r)·(1 + √(1 + (r/2L)²))} + r/(2L) − √(1 + (r/2L)²)]
 *  'long-flush' (L >= 4·width):        Ra = ρ/(2S), S = arithmetic mean of anode length and width
 *  'short-flush' / 'bracelet':         Ra = 0.315·ρ/√A, A = exposed anode surface area
 * Reference: DNV-RP-B401 (anode resistance formulas table).
 * @param {object} o
 * @param {string} o.type see above
 * @param {number} o.resistivity seawater resistivity ρ [Ω·m] (≈ 0.2–0.35 Ω·m; depends on temperature and salinity)
 * @param {number} [o.length] anode length L [m]
 * @param {number} [o.radius] equivalent radius r [m] (r = C/2π for non-cylindrical sections, C = periphery)
 * @param {number} [o.width] anode width [m] (long flush)
 * @param {number} [o.area] exposed surface area A [m²] (short flush / bracelet)
 * @returns {number} resistance [Ω]
 */
export function anodeResistance({ type, resistivity, length, radius, width, area }) {
  positive(resistivity, 'resistivity');
  if (type === 'long-slender-standoff' || type === 'short-slender-standoff') {
    positive(length, 'length');
    positive(radius, 'radius');
    if (type === 'long-slender-standoff') {
      need(length >= 4 * radius, 'length must be at least 4·radius for the long slender formula');
      return (resistivity / (2 * Math.PI * length)) * (Math.log((4 * length) / radius) - 1);
    }
    const x = radius / (2 * length);
    const root = Math.sqrt(1 + x * x);
    return (resistivity / (2 * Math.PI * length)) * (Math.log(((2 * length) / radius) * (1 + root)) + x - root);
  }
  if (type === 'long-flush') {
    positive(length, 'length');
    positive(width, 'width');
    return resistivity / (length + width);
  }
  if (type === 'short-flush' || type === 'bracelet') {
    positive(area, 'area');
    return (0.315 * resistivity) / Math.sqrt(area);
  }
  throw new RangeError(`type: unknown anode type "${type}"`);
}

/**
 * Galvanic anode cathodic protection design.
 *   current demand  Ic = Σ Ac·ic·fc (initial, mean, final)
 *   net anode mass  M = Icm·tf·8760/(u·ε)
 *   anode count     N from mass, then checked for current: N·Ia >= Ic with Ia = (E°c − E°a)/Ra, initially with the
 *                   as-installed dimensions and finally with the depleted anode (stand-off anodes: mass m·(1 − u),
 *                   length reduced by 10 %, radius from the remaining volume around the core; flush and bracelet
 *                   anodes: exposed area unchanged)
 *   capacity        N·ma·ε·u >= Icm·tf·8760
 * Reference: DNV-RP-B401 (2010 and later) calculation procedure and design-parameter tables.
 * Defaults (DNV-RP-B401, Al-Zn-In in seawater): ε = 2000 A·h/kg, closed-circuit anode
 * potential −1.05 V, protective potential −0.80 V (Ag/AgCl/seawater); utilisation 0.90 for long slender
 * stand-off, 0.85 long flush, 0.80 short flush and bracelet. Design current densities depend on depth and climatic
 * region (DNV-RP-B401 design current density tables, roughly 0.06–0.43 A/m²) and are supplied per zone by the caller.
 * @param {object} o
 * @param {{area:number, initialDensity:number, meanDensity:number, finalDensity:number, a?:number, b?:number,
 *   name?:string}[]} o.zones surface area [m²], design current densities [A/m²], coating constants (default bare)
 * @param {number} o.designLife tf [yr]
 * @param {object} o.anode
 * @param {string} o.anode.type anode type, see {@link anodeResistance}
 * @param {number} o.anode.netMass net alloy mass per anode [kg]
 * @param {number} o.anode.length [m]
 * @param {number} [o.anode.radius] initial equivalent radius [m] (stand-off)
 * @param {number} [o.anode.coreRadius=0] radius of the steel core / insert [m] (stand-off)
 * @param {number} [o.anode.width] [m] (long flush)
 * @param {number} [o.anode.area] exposed area [m²] (short flush / bracelet; taken as constant)
 * @param {number} [o.anode.density=2750] alloy density [kg/m³]
 * @param {number} [o.anode.lengthReduction=0.10] fractional length loss of a depleted stand-off anode [-]
 * @param {number} o.resistivity seawater resistivity [Ω·m]
 * @param {number} [o.capacity=2000] electrochemical capacity ε [A·h/kg]
 * @param {number} [o.utilisation] utilisation factor u [-]; default by anode type
 * @param {number} [o.anodePotential=-1.05] closed-circuit anode potential [V]
 * @param {number} [o.protectivePotential=-0.80] design protective potential [V]
 * @returns {{currentInitial:number, currentMean:number, currentFinal:number, requiredMass:number,
 *   anodesForMass:number, anodesForInitialCurrent:number, anodesForFinalCurrent:number, anodes:number,
 *   resistanceInitial:number, resistanceFinal:number, anodeCurrentInitial:number, anodeCurrentFinal:number,
 *   totalNetMass:number, governing:string}} currents [A], mass [kg], resistances [Ω]
 */
export function anodeDesign({ zones, designLife, anode, resistivity, capacity = 2000, utilisation, anodePotential = -1.05, protectivePotential = -0.80 }) {
  need(Array.isArray(zones) && zones.length > 0, 'zones must be a non-empty array');
  positive(designLife, 'designLife');
  positive(capacity, 'capacity');
  need(anode && anode.netMass > 0, 'anode.netMass must be positive');
  const uDefault = { 'long-slender-standoff': 0.90, 'short-slender-standoff': 0.90, 'long-flush': 0.85, 'short-flush': 0.80, bracelet: 0.80 }[anode.type];
  need(uDefault !== undefined, `anode.type: unknown anode type "${anode.type}"`);
  const u = utilisation === undefined ? uDefault : utilisation;
  need(u > 0 && u <= 1, 'utilisation must be in (0, 1]');
  const dE = protectivePotential - anodePotential;
  need(dE > 0, 'protectivePotential must be less negative than anodePotential');
  let Ii = 0;
  let Im = 0;
  let If = 0;
  zones.forEach((z, k) => {
    need(z.area >= 0, `zones[${k}].area must be non-negative`);
    const fc = coatingBreakdown({ a: z.a === undefined ? 1 : z.a, b: z.b === undefined ? 0 : z.b, designLife });
    Ii += z.area * z.initialDensity * fc.initial;
    Im += z.area * z.meanDensity * fc.mean;
    If += z.area * z.finalDensity * fc.final;
  });
  const requiredMass = (Im * designLife * 8760) / (u * capacity);
  const standoff = anode.type.endsWith('standoff');
  const Ri = anodeResistance({ type: anode.type, resistivity, length: anode.length, radius: anode.radius, width: anode.width, area: anode.area });
  let Rf = Ri;
  if (standoff) {
    const rho = anode.density === undefined ? 2750 : anode.density;
    const core = anode.coreRadius || 0;
    const Lf = anode.length * (1 - (anode.lengthReduction === undefined ? 0.10 : anode.lengthReduction));
    const volFinal = (anode.netMass * (1 - u)) / rho;
    const rf = Math.sqrt(volFinal / (Math.PI * Lf) + core * core);
    Rf = anodeResistance({ type: Lf >= 4 * rf ? 'long-slender-standoff' : 'short-slender-standoff', resistivity, length: Lf, radius: rf });
  }
  const Iai = dE / Ri;
  const Iaf = dE / Rf;
  const nMass = Math.ceil(requiredMass / anode.netMass - 1e-9);
  const nInit = Math.ceil(Ii / Iai - 1e-9);
  const nFinal = Math.ceil(If / Iaf - 1e-9);
  const n = Math.max(nMass, nInit, nFinal, 0);
  const governing = n === nMass ? 'mass' : n === nFinal ? 'final-current' : 'initial-current';
  return {
    currentInitial: Ii, currentMean: Im, currentFinal: If, requiredMass, anodesForMass: nMass, anodesForInitialCurrent: nInit,
    anodesForFinalCurrent: nFinal, anodes: n, resistanceInitial: Ri, resistanceFinal: Rf, anodeCurrentInitial: Iai,
    anodeCurrentFinal: Iaf, totalNetMass: n * anode.netMass, governing,
  };
}

// ---------------------------------------------------------------------------------------------------------------
// Marine growth
// ---------------------------------------------------------------------------------------------------------------

/**
 * Default marine-growth thickness profile for the Norwegian continental shelf.
 * Reference: NORSOK N-003 (marine growth): latitude 56–59°N: 100 mm from +2 m to −40 m, 50 mm below −40 m;
 * latitude 59–72°N: 60 mm from +2 m to −40 m, 30 mm below; density 1325 kg/m³. Other regions: supply a profile
 * (e.g. ISO 19901-1 regional annexes).
 * @param {object} o
 * @param {number} o.elevation elevation relative to mean water level, positive up [m]
 * @param {number} [o.latitudeDeg=58] latitude [deg N]
 * @param {{zTop:number, zBottom:number, thickness:number}[]} [o.profile] custom profile overriding the default;
 *   zones ordered from the top, elevations [m], thickness [m]
 * @returns {{thickness:number, density:number}} growth thickness [m] and density [kg/m³]
 */
export function marineGrowth({ elevation, latitudeDeg = 58, profile }) {
  need(Number.isFinite(elevation), 'elevation must be finite');
  const density = 1325;
  if (profile) {
    for (const z of profile) if (elevation <= z.zTop && elevation >= z.zBottom) return { thickness: z.thickness, density };
    return { thickness: 0, density };
  }
  need(latitudeDeg >= 56 && latitudeDeg <= 72, 'latitudeDeg outside the range of the default profile (56–72°N); supply a profile');
  if (elevation > 2) return { thickness: 0, density };
  const south = latitudeDeg < 59;
  if (elevation >= -40) return { thickness: south ? 0.100 : 0.060, density };
  return { thickness: south ? 0.050 : 0.030, density };
}

// ---------------------------------------------------------------------------------------------------------------
// Wear
// ---------------------------------------------------------------------------------------------------------------

/**
 * Archard wear: worn volume V = k·F·s/H (dimensionless coefficient k, hardness H) or V = K·F·s with the
 * dimensional wear coefficient K = k/H [m³/(N·m) = 1/Pa].
 * Reference: Archard (1953) J. Appl. Phys. 24, 981–988. k spans roughly 1e-8 (mild, lubricated) to 1e-3 (severe
 * dry adhesive wear); it must come from tests for the contact pair.
 * @param {object} o
 * @param {number} o.force normal contact force [N]
 * @param {number} o.slidingDistance accumulated sliding distance [m]
 * @param {number} [o.wearCoefficient] dimensional K [1/Pa]
 * @param {number} [o.k] dimensionless coefficient [-] (with hardness)
 * @param {number} [o.hardness] hardness of the softer surface [Pa]
 * @param {number} [o.contactArea] apparent contact area [m²]; when given the mean wear depth is returned
 * @returns {{volume:number, depth:(number|null)}} worn volume [m³], mean depth [m]
 */
export function archardWear({ force, slidingDistance, wearCoefficient, k, hardness, contactArea }) {
  need(force >= 0, 'force must be non-negative');
  need(slidingDistance >= 0, 'slidingDistance must be non-negative');
  let K = wearCoefficient;
  if (K === undefined) {
    positive(k, 'k');
    positive(hardness, 'hardness');
    K = k / hardness;
  }
  need(K >= 0, 'wearCoefficient must be non-negative');
  const volume = K * force * slidingDistance;
  if (contactArea === undefined) return { volume, depth: null };
  positive(contactArea, 'contactArea');
  return { volume, depth: volume / contactArea };
}

/**
 * Interlink wear of mooring chain by Archard's law. Sliding per tension cycle is the relative rotation of
 * adjacent links times the bar radius, s = Δθ·d/2; the worn volume is spread over the interlink contact patch
 * (default area d²·contactFactor) to give a loss of bar diameter at the crown.
 * Reference: Archard (1953). The application to interlink wear is an engineering idealisation: the wear
 * coefficient must be calibrated, and the default contact patch is a geometric assumption, not a code value.
 * @param {object} o
 * @param {number} o.tension mean line tension carried through the contact [N]
 * @param {number} o.diameter bar diameter d [m]
 * @param {number} o.rotationRange relative interlink rotation range per cycle Δθ [rad]
 * @param {number} o.cycles number of cycles [cycles]
 * @param {number} o.wearCoefficient dimensional Archard coefficient K [1/Pa]
 * @param {number} [o.contactFactor=0.25] contact patch area / d² [-]
 * @returns {{volume:number, depth:number, slidingDistance:number, remainingDiameter:number}} [m³], [m], [m], [m]
 */
export function chainInterlinkWear({ tension, diameter, rotationRange, cycles, wearCoefficient, contactFactor = 0.25 }) {
  positive(diameter, 'diameter');
  need(rotationRange >= 0, 'rotationRange must be non-negative');
  need(cycles >= 0, 'cycles must be non-negative');
  positive(contactFactor, 'contactFactor');
  // Each cycle slides over the range and back.
  const s = 2 * rotationRange * (diameter / 2) * cycles;
  const { volume, depth } = archardWear({ force: tension, slidingDistance: s, wearCoefficient, contactArea: contactFactor * diameter * diameter });
  return { volume, depth, slidingDistance: s, remainingDiameter: Math.max(0, diameter - depth) };
}

// ---------------------------------------------------------------------------------------------------------------
// Generic deterioration
// ---------------------------------------------------------------------------------------------------------------

/**
 * Exponential / Weibull-type deterioration of a condition or capacity index: c(t) = c0·exp(−(t/η)^β)
 * (β = 1 is plain exponential decay), with the time at which it reaches a threshold.
 * Reference: generic deterioration curve of Weibull survival form; parameters come from condition data.
 * @param {object} o
 * @param {number} o.time age [yr]
 * @param {number} [o.c0=1] initial condition [-]
 * @param {number} o.eta characteristic time η [yr]
 * @param {number} [o.beta=1] shape β [-]
 * @param {number} [o.threshold] intervention threshold [-], 0 < threshold < c0
 * @returns {{condition:number, timeToThreshold:(number|null), remainingLife:(number|null)}} times [yr]
 */
export function weibullDeterioration({ time, c0 = 1, eta, beta = 1, threshold }) {
  need(time >= 0, 'time must be non-negative');
  positive(c0, 'c0');
  positive(eta, 'eta');
  positive(beta, 'beta');
  const condition = c0 * Math.exp(-((time / eta) ** beta));
  if (threshold === undefined) return { condition, timeToThreshold: null, remainingLife: null };
  need(threshold > 0 && threshold < c0, 'threshold must be in (0, c0)');
  const tt = eta * Math.log(c0 / threshold) ** (1 / beta);
  return { condition, timeToThreshold: tt, remainingLife: Math.max(0, tt - time) };
}

/**
 * Condition-state Markov deterioration: discrete states (1 = as new … n = failed/worst) with an annual
 * transition matrix P. Returns the state distribution over time, the expected state and the expected time to
 * first reach the worst state from each state.
 * Reference: Markov deterioration models of bridge and asset management systems (e.g. Jiang & Sinha 1989;
 * Frangopol, Kallen & van Noortwijk 2004).
 * @param {object} o
 * @param {number[][]} o.P annual transition matrix, rows sum to 1, last state absorbing for the life estimate
 * @param {number[]} o.initial initial state distribution (sums to 1)
 * @param {number} o.years projection horizon [yr], integer
 * @param {number} [o.failureProbabilityLimit] optional limit on P[worst state]; the first year exceeding it is returned
 * @returns {{distribution:number[][], expectedState:number[], pWorst:number[], meanTimeToWorst:number[],
 *   yearLimitExceeded:(number|null)}} distribution[y][state]; expectedState 1-based; meanTimeToWorst [yr] per
 *   starting state (0 for the worst state, Infinity if it cannot be reached)
 */
export function markovDeterioration({ P, initial, years, failureProbabilityLimit }) {
  need(Array.isArray(P) && P.length >= 2, 'P must be a square matrix of size >= 2');
  const n = P.length;
  P.forEach((row, i) => {
    need(row.length === n, 'P must be square');
    let s = 0;
    for (const v of row) { need(v >= 0, 'P entries must be non-negative'); s += v; }
    need(Math.abs(s - 1) < 1e-9, `P row ${i} must sum to 1`);
  });
  need(Array.isArray(initial) && initial.length === n, 'initial must have one entry per state');
  need(Number.isInteger(years) && years >= 0, 'years must be a non-negative integer');
  let p = initial.slice();
  const distribution = [p.slice()];
  for (let y = 0; y < years; y++) {
    const q = new Array(n).fill(0);
    for (let i = 0; i < n; i++) if (p[i] !== 0) for (let j = 0; j < n; j++) q[j] += p[i] * P[i][j];
    p = q;
    distribution.push(p.slice());
  }
  const expectedState = distribution.map((d) => d.reduce((s, v, i) => s + v * (i + 1), 0));
  const pWorst = distribution.map((d) => d[n - 1]);
  // Mean first-passage time to the worst state: (I − P_TT)·τ = 1 on the other states, by Gauss–Seidel
  // (P_TT is sub-stochastic, so the iteration converges when the worst state is reachable).
  const tau = new Array(n).fill(0);
  let converged = false;
  for (let it = 0; it < 200000 && !converged; it++) {
    let diff = 0;
    for (let i = 0; i < n - 1; i++) {
      let s = 1;
      for (let j = 0; j < n - 1; j++) if (j !== i) s += P[i][j] * tau[j];
      const den = 1 - P[i][i];
      const v = den > 1e-14 ? s / den : Infinity;
      diff = Math.max(diff, Math.abs(v - tau[i]) / Math.max(1, Math.abs(v)));
      tau[i] = v;
    }
    converged = diff < 1e-12 || !Number.isFinite(diff);
  }
  const meanTimeToWorst = tau.map((v, i) => (i === n - 1 ? 0 : converged && Number.isFinite(v) ? v : Infinity));
  let yearLimitExceeded = null;
  if (failureProbabilityLimit !== undefined) {
    const k = pWorst.findIndex((v) => v > failureProbabilityLimit);
    yearLimitExceeded = k >= 0 ? k : null;
  }
  return { distribution, expectedState, pWorst, meanTimeToWorst, yearLimitExceeded };
}
