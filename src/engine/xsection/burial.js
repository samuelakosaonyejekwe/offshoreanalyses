/**
 * Protection of cables, umbilicals and small pipelines by burial and rock cover: burial protection index,
 * probabilistic anchor-strike assessment in the structure of the Carbon Trust Cable Burial Risk Assessment
 * (CBRA) method, fishing-gear penetration, trenching tool performance, rock berms, mattresses and sandwave
 * exposure.
 *
 * Every empirical factor is an argument; the tables shipped here are flagged generic and editable.
 */

function pos(v, name) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}
function fin(v, name) {
  if (!(typeof v === 'number' && Number.isFinite(v))) throw new RangeError(`${name} must be a finite number`);
  return v;
}
function nonneg(v, name) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v >= 0)) throw new RangeError(`${name} must be a non-negative finite number`);
  return v;
}

/**
 * Burial Protection Index classes (Mole, Featherstone & Winter 1997; Allan 1998): the BPI relates burial depth
 * to the depth needed in the same soil for a given level of protection.
 *   BPI 1  protection against normal fishing gear only (areas where anchoring is unlikely)
 *   BPI 2  protection against anchors up to about 2 tonnes (vessels other than larger ships)
 *   BPI 3  protection against anchors of all but the largest ships
 */
export const BPI_CLASSES = Object.freeze([
  Object.freeze({ bpi: 1, threat: 'fishing', description: 'Normal fishing gear only' }),
  Object.freeze({ bpi: 2, threat: 'anchor-small', description: 'Anchors up to approximately 2 tonnes' }),
  Object.freeze({ bpi: 3, threat: 'anchor-large', description: 'Anchors of all but the largest ships' }),
]);

/**
 * Burial Protection Index of a burial depth: BPI = burial depth / depth giving BPI = 1 in the same soil, and the
 * depth required for a target BPI. The reference depth for BPI = 1 depends on soil type and strength (it is
 * larger in soft clay than in stiff clay or dense sand) and must be read from the published BPI chart or a
 * project-specific study: it is a required argument.
 * Reference: Mole, Featherstone & Winter (1997) "Cable protection - solutions through new installation and
 * burial approaches", SubOptic '97; Allan (1998) "Selecting appropriate cable burial depths - a methodology",
 * IBC Submarine Communications conference.
 * Validity: BPI is an index, linear in depth by definition; it is superseded for wind-farm cables by
 * risk-based methods ({@link cbraStrikeProbability}).
 * @param {object} o
 * @param {number} o.referenceDepth burial depth giving BPI = 1 in this soil, m
 * @param {number} [o.burialDepth] achieved depth, m
 * @param {number} [o.targetBpi] required index
 * @returns {{bpi:number|null, requiredDepth:number|null, protectedAgainst:string|null}} -, m
 */
export function burialProtectionIndex({ referenceDepth, burialDepth, targetBpi }) {
  pos(referenceDepth, 'referenceDepth');
  const bpi = burialDepth === undefined ? null : nonneg(burialDepth, 'burialDepth') / referenceDepth;
  let protectedAgainst = null;
  if (bpi !== null) protectedAgainst = bpi >= 3 ? BPI_CLASSES[2].description : bpi >= 2 ? BPI_CLASSES[1].description : bpi >= 1 ? BPI_CLASSES[0].description : 'Below BPI 1';
  return { bpi, requiredDepth: targetBpi === undefined ? null : pos(targetBpi, 'targetBpi') * referenceDepth, protectedAgainst };
}

/**
 * Anchor penetration depth for one vessel/anchor class.
 *     penetration = flukeLength x sin(flukeAngle) x soilFactor
 * The fluke-tip depth of a fully tripped stockless anchor on a firm seabed is flukeLength x sin(flukeAngle)
 * (geometry); `soilFactor` scales it for the soil (close to 1 in sand and stiff clay, greater than 1 in soft clay
 * where the whole anchor sinks; the CBRA guidance tabulates values by soil type - argument). When the fluke
 * length is not known it can be estimated by geometric similarity, flukeLength = flukeCoefficient x
 * (anchorMass)^(1/3), with the coefficient fitted to the anchor type's catalogue dimensions (argument, m/kg^(1/3)).
 * Reference: Carbon Trust (2015) Cable Burial Risk Assessment Methodology (anchor penetration); geometry.
 * @param {object} o
 * @param {number} [o.flukeLength] m
 * @param {number} [o.anchorMass] kg (with flukeCoefficient)
 * @param {number} [o.flukeCoefficient] m/kg^(1/3)
 * @param {number} [o.flukeAngleDeg=45] fluke opening angle, deg (stockless anchors open to about 35-45 deg)
 * @param {number} [o.soilFactor=1] -
 * @returns {{flukeLength:number, penetration:number}} m
 */
export function anchorPenetration({ flukeLength, anchorMass, flukeCoefficient, flukeAngleDeg = 45, soilFactor = 1 }) {
  let Lf = flukeLength;
  if (Lf === undefined) Lf = pos(flukeCoefficient, 'flukeCoefficient') * Math.cbrt(pos(anchorMass, 'anchorMass'));
  pos(Lf, 'flukeLength'); pos(soilFactor, 'soilFactor');
  if (!(flukeAngleDeg > 0 && flukeAngleDeg <= 90)) throw new RangeError('flukeAngleDeg must be in (0, 90]');
  return { flukeLength: Lf, penetration: Lf * Math.sin(flukeAngleDeg * Math.PI / 180) * soilFactor };
}

/**
 * Annual probability of anchor strike on a cable route section, in the structure of the CBRA method:
 *     P_strike = P_traffic x P_wd x sum over vessels [ D_ship / (V_ship x 8760 h) x P_incident ]
 * where for each vessel crossing the section D_ship is the distance travelled within the zone from which a
 * dragged or dropped anchor can reach the cable, V_ship its speed (so D/V is the exposure time), P_incident the
 * probability per vessel-year at sea of an incident leading to anchor deployment (e.g. loss of propulsion),
 * P_wd a modifier for water depth / seabed (anchoring feasibility) and P_traffic a traffic modifier. Here the
 * vessels are grouped in classes with `crossings` per year (from AIS data); only classes whose anchor
 * penetration exceeds the burial depth contribute. Return period = 1 / P_strike.
 * Reference: Carbon Trust (2015) Cable Burial Risk Assessment Methodology - Guidance for the Preparation of
 * Cable Burial Depth of Lowering Specification (probabilistic assessment). All modifiers are arguments.
 * Validity: rare-event sum (P << 1); independent crossings.
 * @param {object} o
 * @param {{name?:string, crossings:number, speed:number, hazardDistance:number, pIncident:number,
 *   penetration:number, pWd?:number}[]} o.classes crossings per year; speed m/s; hazardDistance D_ship m;
 *   pIncident per vessel-year; anchor penetration m; class-specific depth modifier (default 1)
 * @param {number} [o.burialDepth=0] depth of lowering to top of cable, m
 * @param {number} [o.pTraffic=1] traffic modifier
 * @param {number} [o.pWd=1] water-depth modifier applied to all classes
 * @returns {{pStrike:number, returnPeriod:number, classes:{name:string, exposureHours:number,
 *   contribution:number, threatens:boolean}[]}} per year, years
 */
export function cbraStrikeProbability({ classes, burialDepth = 0, pTraffic = 1, pWd = 1 }) {
  if (!Array.isArray(classes)) throw new RangeError('classes must be an array');
  nonneg(burialDepth, 'burialDepth'); nonneg(pTraffic, 'pTraffic'); nonneg(pWd, 'pWd');
  let p = 0;
  const rows = classes.map((c, i) => {
    nonneg(c.crossings, `classes[${i}].crossings`); pos(c.speed, `classes[${i}].speed`); nonneg(c.hazardDistance, `classes[${i}].hazardDistance`);
    nonneg(c.pIncident, `classes[${i}].pIncident`); nonneg(c.penetration, `classes[${i}].penetration`);
    const hours = c.crossings * c.hazardDistance / c.speed / 3600;
    const threatens = c.penetration > burialDepth;
    const contribution = threatens ? pTraffic * pWd * (c.pWd ?? 1) * hours / 8760 * c.pIncident : 0;
    p += contribution;
    return { name: c.name ?? `class ${i + 1}`, exposureHours: hours, contribution, threatens };
  });
  return { pStrike: p, returnPeriod: p > 0 ? 1 / p : Infinity, classes: rows };
}

/**
 * Depth-of-lowering recommendation: the smallest burial depth for which the CBRA strike probability does not
 * exceed the target (the candidate depths are the anchor penetrations of the classes, since the probability only
 * changes there), plus a margin; never less than the fishing-gear penetration plus margin.
 * Reference: Carbon Trust (2015) CBRA (depth of lowering specification).
 * @param {object} o
 * @param {object[]} o.classes as for {@link cbraStrikeProbability}
 * @param {number} o.targetProbability acceptable annual strike probability (e.g. 1e-4 to 1e-5 per section: project
 *   risk criterion)
 * @param {number} [o.fishingPenetration=0] m
 * @param {number} [o.margin=0] safety margin added to the governing penetration, m
 * @param {number} [o.pTraffic=1]
 * @param {number} [o.pWd=1]
 * @returns {{recommendedDepth:number, governing:'anchor'|'fishing'|'none', pStrike:number, returnPeriod:number,
 *   pStrikeUnburied:number}} m, -, 1/yr, yr, 1/yr
 */
export function depthOfLowering({ classes, targetProbability, fishingPenetration = 0, margin = 0, pTraffic = 1, pWd = 1 }) {
  pos(targetProbability, 'targetProbability'); nonneg(fishingPenetration, 'fishingPenetration'); nonneg(margin, 'margin');
  const cands = [0, ...classes.map((c) => c.penetration)].sort((a, b) => a - b);
  let depth = cands[cands.length - 1];
  for (const d of cands) {
    if (cbraStrikeProbability({ classes, burialDepth: d, pTraffic, pWd }).pStrike <= targetProbability) { depth = d; break; }
  }
  const gov = depth > fishingPenetration ? 'anchor' : fishingPenetration > 0 ? 'fishing' : 'none';
  const base = Math.max(depth, fishingPenetration);
  const res = cbraStrikeProbability({ classes, burialDepth: base, pTraffic, pWd });
  return {
    recommendedDepth: base > 0 ? base + margin : 0, governing: gov, pStrike: res.pStrike, returnPeriod: res.returnPeriod,
    pStrikeUnburied: cbraStrikeProbability({ classes, burialDepth: 0, pTraffic, pWd }).pStrike,
  };
}

/**
 * GENERIC, editable table of typical fishing-gear seabed penetration depths (m) by gear type and soil. The
 * values are rounded indicative figures of the order reported in fishing-gear interaction studies (penetration
 * of trawl doors, beam-trawl shoes/tickler chains and dredge teeth is generally below 0.3 m and larger in soft
 * mud than in sand); they are not taken from a single source and must be confirmed for the fishery concerned.
 */
export const GENERIC_FISHING_PENETRATION = Object.freeze({
  generic: true,
  otterTrawl: Object.freeze({ sand: 0.2, mud: 0.3, stiffClay: 0.1, gravel: 0.1 }),
  beamTrawl: Object.freeze({ sand: 0.1, mud: 0.3, stiffClay: 0.1, gravel: 0.1 }),
  scallopDredge: Object.freeze({ sand: 0.15, mud: 0.2, stiffClay: 0.1, gravel: 0.15 }),
  hydraulicDredge: Object.freeze({ sand: 0.3, mud: 0.3, stiffClay: 0.2, gravel: 0.2 }),
});

/**
 * Fishing-gear penetration lookup with a safety factor.
 * @param {object} o
 * @param {string} o.gear key of the table
 * @param {string} o.soil key of the gear row
 * @param {number} [o.safetyFactor=1]
 * @param {object} [o.table=GENERIC_FISHING_PENETRATION] editable table { gear: { soil: depth [m] } }
 * @returns {{penetration:number, design:number, generic:boolean}} m
 */
export function fishingGearPenetration({ gear, soil, safetyFactor = 1, table = GENERIC_FISHING_PENETRATION }) {
  pos(safetyFactor, 'safetyFactor');
  const row = table[gear];
  if (!row || typeof row !== 'object') throw new RangeError(`gear '${gear}' is not in the table`);
  const d = row[soil];
  if (!(typeof d === 'number' && d >= 0)) throw new RangeError(`soil '${soil}' is not in the table for gear '${gear}'`);
  return { penetration: d, design: d * safetyFactor, generic: table.generic === true };
}

/**
 * Progress rate of a trenching tool from an energy balance:
 *     v = efficiency x power / (specificEnergy x trenchArea)
 * with the trench cross-section A = width x depth (or the supplied area) and the specific energy e (J per m^3 of
 * soil disaggregated/fluidised/cut). For water jetting in clay e scales with the undrained shear strength
 * (e = jetFactor x su, the factor being a tool characteristic from trials); for mechanical chain cutters in
 * weak rock e scales with the unconfined compressive strength. Both scalings are available through
 * `strength` x `strengthFactor`.
 * Reference: specific-energy concept of excavation (Teale 1965, Int. J. Rock Mech. Min. Sci. 2, 57-73) applied
 * as a documented energy model; not a tool-specific performance curve.
 * Validity: steady progress, power-limited (not tow-force- or stability-limited).
 * @param {object} o
 * @param {number} o.power installed jetting or cutting power, W
 * @param {number} [o.efficiency=1] fraction of the power reaching the soil
 * @param {number} [o.specificEnergy] e, J/m^3
 * @param {number} [o.strength] su or UCS, Pa (used when specificEnergy is omitted)
 * @param {number} [o.strengthFactor] e / strength, -
 * @param {number} [o.trenchWidth] m
 * @param {number} [o.trenchDepth] m
 * @param {number} [o.trenchArea] m^2 (instead of width x depth)
 * @returns {{rate:number, ratePerHour:number, specificEnergy:number, trenchArea:number}} m/s, m/h, J/m^3, m^2
 */
export function trenchingRate({ power, efficiency = 1, specificEnergy, strength, strengthFactor, trenchWidth, trenchDepth, trenchArea }) {
  pos(power, 'power');
  if (!(efficiency > 0 && efficiency <= 1)) throw new RangeError('efficiency must be in (0, 1]');
  const e = specificEnergy !== undefined ? pos(specificEnergy, 'specificEnergy') : pos(strength, 'strength') * pos(strengthFactor, 'strengthFactor');
  const A = trenchArea !== undefined ? pos(trenchArea, 'trenchArea') : pos(trenchWidth, 'trenchWidth') * pos(trenchDepth, 'trenchDepth');
  const rate = efficiency * power / (e * A);
  return { rate, ratePerHour: rate * 3600, specificEnergy: e, trenchArea: A };
}

/**
 * Plough tow force.
 *  - Empirical form of Cathie & Wintgens (2001, OTC 13145):
 *        sand:  F = Cw W' + Cs gamma' D^3 + Cd v D^2
 *        clay:  F = Cw W' + Cc su D^2
 *    (W' plough submerged weight, D trench depth, v speed, gamma' submerged unit weight, su undrained shear
 *    strength). The coefficients Cw (interface friction), Cs (static, depends on relative density), Cd (dynamic
 *    rate term from dilation and permeability, N s/m^3) and Cc are plough- and soil-specific: required arguments.
 *  - First-principles static estimate reported alongside: Rankine passive wedge on the share front of width b,
 *        sand:  P = 0.5 Kp gamma' D^2 b,  Kp = tan^2(45 + phi/2)
 *        clay:  P = (0.5 gamma' D^2 + 2 su D) b
 *    plus base friction mu W'.
 * Reference: Cathie & Wintgens (2001) "Pipeline trenching using plows: performance and geotechnical hazards",
 * OTC 13145 (functional form); Rankine (1857) passive pressure.
 * Validity: steady ploughing in uniform soil; the Rankine estimate ignores side-wall shear and three-dimensional
 * wedge effects (lower bound of the static part).
 * @param {object} o
 * @param {'sand'|'clay'} o.soil
 * @param {number} o.depth D, m
 * @param {number} o.ploughWeight W' submerged, N
 * @param {number} o.gammaSub gamma', N/m^3
 * @param {number} [o.speed=0] v, m/s
 * @param {number} [o.Cw] -
 * @param {number} [o.Cs] - (sand)
 * @param {number} [o.Cd] N s/m^3 (sand)
 * @param {number} [o.Cc] - (clay)
 * @param {number} [o.su] Pa (clay)
 * @param {number} [o.phiDeg] friction angle, deg (sand, Rankine)
 * @param {number} [o.shareWidth] b, m (Rankine)
 * @param {number} [o.mu] base friction coefficient (Rankine; default Cw if given)
 * @returns {{empirical:number|null, staticPart:number|null, ratePart:number|null, rankine:number|null}} N
 */
export function ploughTowForce({ soil, depth, ploughWeight, gammaSub, speed = 0, Cw, Cs, Cd, Cc, su, phiDeg, shareWidth, mu }) {
  if (soil !== 'sand' && soil !== 'clay') throw new RangeError("soil must be 'sand' or 'clay'");
  const D = pos(depth, 'depth'); nonneg(ploughWeight, 'ploughWeight'); pos(gammaSub, 'gammaSub'); nonneg(speed, 'speed');
  let empirical = null, staticPart = null, ratePart = null, rankine = null;
  if (soil === 'sand') {
    if (Cw !== undefined && Cs !== undefined) {
      staticPart = nonneg(Cw, 'Cw') * ploughWeight + nonneg(Cs, 'Cs') * gammaSub * D ** 3;
      ratePart = Cd !== undefined ? nonneg(Cd, 'Cd') * speed * D * D : 0;
      empirical = staticPart + ratePart;
    }
    if (phiDeg !== undefined && shareWidth !== undefined) {
      if (!(phiDeg > 0 && phiDeg < 90)) throw new RangeError('phiDeg must be in (0, 90)');
      const Kp = Math.tan((45 + phiDeg / 2) * Math.PI / 180) ** 2;
      rankine = 0.5 * Kp * gammaSub * D * D * pos(shareWidth, 'shareWidth') + (mu ?? Cw ?? 0) * ploughWeight;
    }
  } else {
    pos(su, 'su');
    if (Cw !== undefined && Cc !== undefined) {
      staticPart = nonneg(Cw, 'Cw') * ploughWeight + nonneg(Cc, 'Cc') * su * D * D; ratePart = 0; empirical = staticPart;
    }
    if (shareWidth !== undefined) rankine = (0.5 * gammaSub * D * D + 2 * su * D) * pos(shareWidth, 'shareWidth') + (mu ?? Cw ?? 0) * ploughWeight;
  }
  if (empirical === null && rankine === null) throw new RangeError('give the empirical coefficients (Cw with Cs or Cc) and/or the Rankine inputs (shareWidth, phiDeg)');
  return { empirical, staticPart, ratePart, rankine };
}

/**
 * Time for an open trench to fill naturally: t = trenchArea / (transportRate x trapEfficiency), with the
 * sediment transport rate per unit width normal to the trench (m^2/s, i.e. m^3 per metre per second, from a
 * sediment-transport calculation) and the fraction trapped. Mechanical backfill (plough backfill pass, jetting
 * re-settlement) is immediate but typically leaves loose, lower-strength, higher-thermal-resistivity cover.
 * Reference: sediment continuity (Exner equation integrated over the trench).
 * @param {object} o
 * @param {number} o.trenchArea open cross-section, m^2
 * @param {number} o.transportRate m^2/s
 * @param {number} [o.trapEfficiency=1] (0, 1]
 * @returns {{time:number, timeDays:number}} s, days
 */
export function naturalInfillTime({ trenchArea, transportRate, trapEfficiency = 1 }) {
  pos(trenchArea, 'trenchArea'); pos(transportRate, 'transportRate');
  if (!(trapEfficiency > 0 && trapEfficiency <= 1)) throw new RangeError('trapEfficiency must be in (0, 1]');
  const t = trenchArea / (transportRate * trapEfficiency);
  return { time: t, timeDays: t / 86400 };
}

/**
 * Rock berm over a cable or pipe: trapezoidal geometry and tonnage take-off.
 * Height h = product diameter + cover above the product (plus any under-layer); cross-section
 * A = h (crestWidth + slope h) for side slopes 1:slope (vertical:horizontal); base width = crest + 2 slope h; rock
 * volume per metre = A - product area (if not in a trench); mass = volume x bulk density x (1 + overdump).
 * The cover must also satisfy any impact/penetration absorption thickness (anchor or trawl board; argument).
 * Reference: geometry; rock berm design practice per CIRIA C683 The Rock Manual (2007) in intent.
 * @param {object} o
 * @param {number} o.productDiameter m
 * @param {number} o.cover rock cover above top of product, m
 * @param {number} o.crestWidth m
 * @param {number} [o.slope=3] horizontal per unit vertical (1:3 typical for fall-pipe berms; argument)
 * @param {number} o.length berm length, m
 * @param {number} [o.bulkDensity] kg/m^3; or
 * @param {number} [o.rockDensity=2650] solid density, kg/m^3
 * @param {number} [o.porosity=0.4] -
 * @param {number} [o.overdump=0] fractional allowance for placement tolerance and losses
 * @param {number} [o.requiredAbsorptionCover=0] minimum cover for impact energy absorption, m
 * @returns {{height:number, baseWidth:number, area:number, volume:number, mass:number, tonnes:number,
 *   coverUtilisation:number}} m, m, m^2, m^3, kg, t, required/provided cover
 */
export function rockBerm({ productDiameter, cover, crestWidth, slope = 3, length, bulkDensity, rockDensity = 2650, porosity = 0.4, overdump = 0, requiredAbsorptionCover = 0 }) {
  nonneg(productDiameter, 'productDiameter'); pos(cover, 'cover'); nonneg(crestWidth, 'crestWidth'); pos(slope, 'slope'); pos(length, 'length');
  nonneg(overdump, 'overdump'); nonneg(requiredAbsorptionCover, 'requiredAbsorptionCover');
  if (!(porosity >= 0 && porosity < 1)) throw new RangeError('porosity must be in [0, 1)');
  const h = productDiameter + cover;
  const gross = h * (crestWidth + slope * h);
  const area = gross - Math.PI * productDiameter ** 2 / 4;
  const rhoB = bulkDensity !== undefined ? pos(bulkDensity, 'bulkDensity') : pos(rockDensity, 'rockDensity') * (1 - porosity);
  const volume = area * length * (1 + overdump);
  return { height: h, baseWidth: crestWidth + 2 * slope * h, area, volume, mass: volume * rhoB, tonnes: volume * rhoB / 1000, coverUtilisation: requiredAbsorptionCover / cover };
}

/**
 * Rock size for hydraulic stability from the Shields criterion (hook): the bed shear stress tau must not exceed
 * the critical value, tau_cr = psi_cr (rho_s - rho_w) g D50, hence D50 >= tau / (psi_cr (rho_s - rho_w) g).
 * tau may be given directly or as rho_w u*^2 from a shear velocity. The critical Shields parameter is an
 * argument (about 0.03-0.035 for no movement and 0.05-0.055 for limited movement of coarse rock, CIRIA C683
 * - confirm for the design approach).
 * Reference: Shields (1936); CIRIA C683 The Rock Manual (2007).
 * @param {object} o
 * @param {number} [o.bedShear] tau, Pa
 * @param {number} [o.shearVelocity] u*, m/s
 * @param {number} o.shieldsCritical psi_cr, -
 * @param {number} [o.rockDensity=2650] kg/m^3
 * @param {number} [o.rhoW=1025] kg/m^3
 * @param {number} [o.g=9.80665] m/s^2
 * @param {number} [o.d50] provided rock size, m (for the utilisation)
 * @returns {{requiredD50:number, bedShear:number, utilisation:number|null}} m, Pa, -
 */
export function rockSizeShields({ bedShear, shearVelocity, shieldsCritical, rockDensity = 2650, rhoW = 1025, g = 9.80665, d50 }) {
  pos(shieldsCritical, 'shieldsCritical'); pos(rhoW, 'rhoW'); pos(g, 'g');
  if (!(rockDensity > rhoW)) throw new RangeError('rockDensity must exceed rhoW');
  const tau = bedShear !== undefined ? nonneg(bedShear, 'bedShear') : rhoW * nonneg(shearVelocity, 'shearVelocity') ** 2;
  const req = tau / (shieldsCritical * (rockDensity - rhoW) * g);
  return { requiredD50: req, bedShear: tau, utilisation: d50 === undefined ? null : req / pos(d50, 'd50') };
}

/**
 * Number of concrete mattresses to cover a length: n = ceil((length + endOverrun) / (mattressLength - overlap))
 * per row, times rows = ceil(coverWidth / mattressWidth); total submerged weight for vessel planning.
 * @param {object} o
 * @param {number} o.length protected length, m
 * @param {number} o.mattressLength dimension along the product, m
 * @param {number} [o.overlap=0] overlap (positive) or gap (negative) between adjacent mattresses, m
 * @param {number} [o.endOverrun=0] extra length beyond the ends, m
 * @param {number} [o.coverWidth] required width, m
 * @param {number} [o.mattressWidth] m
 * @param {number} [o.mattressMass] kg each
 * @returns {{perRow:number, rows:number, count:number, totalMass:number|null}} -, -, -, kg
 */
export function mattressCount({ length, mattressLength, overlap = 0, endOverrun = 0, coverWidth, mattressWidth, mattressMass }) {
  pos(length, 'length'); pos(mattressLength, 'mattressLength'); fin(overlap, 'overlap'); nonneg(endOverrun, 'endOverrun');
  if (!(overlap < mattressLength)) throw new RangeError('overlap must be smaller than mattressLength');
  const perRow = Math.ceil((length + endOverrun) / (mattressLength - overlap) - 1e-9);
  const rows = coverWidth !== undefined ? Math.ceil(pos(coverWidth, 'coverWidth') / pos(mattressWidth, 'mattressWidth') - 1e-9) : 1;
  const count = perRow * rows;
  return { perRow, rows, count, totalMass: mattressMass === undefined ? null : count * pos(mattressMass, 'mattressMass') };
}

/**
 * Thermal record of a burial condition for the electrical (ampacity) engine: burial depth to the cable axis and
 * the thermal resistivity of the backfill and native soil. Also returns the external thermal resistance of a
 * single isolated buried cable, T4 = (rho_T / 2 pi) ln(u + sqrt(u^2 - 1)), u = 2 L / De (IEC 60287-2-1, single isolated buried cable),
 * so the effect of the depth of lowering on rating can be seen directly.
 * @param {object} o
 * @param {number} o.depthToTop depth of lowering to top of cable, m
 * @param {number} o.outerDiameter De, m
 * @param {number} o.soilResistivity native soil thermal resistivity, K m/W
 * @param {number} [o.backfillResistivity] K m/W (default: soil)
 * @returns {{burialDepthAxis:number, soilResistivity:number, backfillResistivity:number, u:number,
 *   T4:number}} m, K m/W, K m/W, -, K m/W (with the backfill resistivity)
 */
export function burialThermalRecord({ depthToTop, outerDiameter, soilResistivity, backfillResistivity }) {
  nonneg(depthToTop, 'depthToTop'); pos(outerDiameter, 'outerDiameter'); pos(soilResistivity, 'soilResistivity');
  const rb = backfillResistivity === undefined ? soilResistivity : pos(backfillResistivity, 'backfillResistivity');
  const L = depthToTop + outerDiameter / 2;
  const u = 2 * L / outerDiameter;
  return { burialDepthAxis: L, soilResistivity, backfillResistivity: rb, u, T4: rb / (2 * Math.PI) * Math.log(u + Math.sqrt(u * u - 1)) };
}

/**
 * Sandwave mobility exposure (geometric model). The seabed is a sinusoidal bedform of height H (trough to crest)
 * and wavelength lambda migrating at celerity c without change of shape:
 *     z(x, t) = (H/2) cos(2 pi (x - c t)/lambda).
 * A cable installed at time 0 with cover d0 at a point whose phase is phi0 (0 = crest, pi = trough) has cover
 *     d(t) = d0 + (H/2) [cos(phi0 - 2 pi c t/lambda) - cos(phi0)].
 * The first time the cover falls to the target d_t is found in closed form; if the cover never falls that far
 * the time is Infinity. The worst case is a cable laid under a crest (phi0 = 0):
 * t = lambda/(2 pi c) acos(1 - 2 (d0 - d_t)/H). Minimum cover over a full cycle is d0 - (H/2)(1 + cos(phi0)).
 * Reference: kinematic bedform migration (geometry); migration rates and heights from repeat bathymetric surveys.
 * Validity: steady migration of an unchanged profile; ignores megaripple superposition and bedform growth.
 * @param {object} o
 * @param {number} o.height H, m
 * @param {number} o.wavelength lambda, m
 * @param {number} o.migrationRate c, m/yr (or any time unit; the result uses the same unit)
 * @param {number} o.initialCover d0, m
 * @param {number} [o.targetCover=0] d_t, m
 * @param {number} [o.phase=0] phi0, rad
 * @returns {{time:number, minimumCover:number, exposed:boolean, period:number}} time until the cover reaches the
 *   target (same time unit as migrationRate), m, whether the target is ever reached, migration period lambda/c
 */
export function sandwaveExposure({ height, wavelength, migrationRate, initialCover, targetCover = 0, phase = 0 }) {
  const H = pos(height, 'height'), lam = pos(wavelength, 'wavelength'), c = pos(migrationRate, 'migrationRate');
  nonneg(initialCover, 'initialCover'); fin(targetCover, 'targetCover'); fin(phase, 'phase');
  const period = lam / c;
  const minimumCover = initialCover - H / 2 * (1 + Math.cos(phase));
  if (initialCover <= targetCover) return { time: 0, minimumCover, exposed: true, period };
  const C = Math.cos(phase) - 2 * (initialCover - targetCover) / H;
  if (C < -1) return { time: Infinity, minimumCover, exposed: false, period };
  // theta(t) = phase - 2 pi c t/lambda decreases; first theta below phase with cos(theta) = C
  const a = Math.acos(Math.max(C, -1));
  const p0 = Math.atan2(Math.sin(phase), Math.cos(phase)); // (-pi, pi]
  let best = Infinity;
  for (let k = -2; k <= 1; k++) {
    for (const th of [a + 2 * Math.PI * k, -a + 2 * Math.PI * k]) {
      const d = p0 - th;
      if (d > 1e-14 && d < best) best = d;
    }
  }
  return { time: best * lam / (2 * Math.PI * c), minimumCover, exposed: true, period };
}
