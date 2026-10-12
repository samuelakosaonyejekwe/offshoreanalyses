/**
 * Pipeline, cable and umbilical installation statics: S-lay over a circular stinger, J-lay, reel-lay bending,
 * cable lay over a chute or lay wheel, abandonment and recovery, shore pull-in, plough tow-force hook,
 * quasi-dynamic tension bounds, the boundary-layer stiffened catenary and DNV-ST-F101 simplified laying criteria.
 *
 * Coordinates for the rod-based analyses: seabed at z = 0, x positive from the vessel towards the laid pipe.
 * All analyses are in-plane and static unless stated.
 */

import { solveRodStatic } from './rod.js';
import { solveScr, catenaryRiser } from './riser.js';
import { brent } from '../math/roots.js';

const DEG = Math.PI / 180;

function num(v, name, { min = -Infinity, max = Infinity, strictMin = false } = {}) {
  if (typeof v !== 'number' || !Number.isFinite(v) || v < min || v > max || (strictMin && v <= min)) {
    throw new RangeError(`pipelay: argument ${name} = ${v} is outside its valid range`);
  }
  return v;
}

/* ------------------------------------------------------------------------------------------------------- */
/* S-lay                                                                                                     */

/* Validated inputs → rod model, mesh bookkeeping and the natural lift-off geometry shared by sLay and sLayRodModel. */
function sLaySetup(o) {
  const { pipe, D, hingeHeight: zH, stingerRadius: R, stingerLength: Ls, tension: T, nElem = 260, rollerStiffness = 1e8 } = o ?? {};
  if (!pipe) throw new RangeError('pipelay: argument pipe is required');
  num(D, 'D', { min: 0, strictMin: true }); num(zH, 'hingeHeight', { min: 0, strictMin: true });
  num(R, 'stingerRadius', { min: 0, strictMin: true }); num(Ls, 'stingerLength', { min: 0, strictMin: true, max: 0.5 * Math.PI * R });
  const w = pipe.w;
  num(T, 'tension', { min: w * zH, strictMin: true });
  const H = T - w * zH, a = H / w;
  // natural lift-off: stinger slope equals the catenary slope at the same elevation
  const thTip = Ls / R;
  const g = (phi) => Math.cos(phi) - H / (H + w * (zH - R * (1 - Math.cos(phi))));
  let phiL = thTip;
  if (g(1e-9) > 0 && g(thTip) < 0) phiL = brent(g, 1e-9, thTip, { tol: 1e-12 });
  const sLift = R * phiL, zL = zH - R * (1 - Math.cos(phiL));
  const sig = Math.sqrt(zL * (zL + 2 * a));
  const Lg = num(o.groundedLength ?? 0.35 * sig, 'groundedLength', { min: 0, strictMin: true });
  const lob = Math.min(1.3 * Ls, sLift + 0.5 * sig);
  const n1 = Math.max(30, Math.round(0.27 * nElem)), n3 = Math.max(10, Math.round(0.13 * nElem)), n2 = Math.max(40, nElem - n1 - n3);
  const sections = [
    { ...pipe, length: lob, nElem: n1 },
    { ...pipe, length: sLift + sig - lob + 0.4 * Lg, nElem: n2 },
    { ...pipe, length: 0.6 * Lg, nElem: n3 },
  ];
  // Starting shape, marched from the tensioner along the stretched element lengths so that it already carries the
  // lay tension Te = H + w·z (frictionless stinger and catenary alike); a stress-free start has no geometric
  // stiffness against the tensioner force and makes the first equilibrium iteration ill-conditioned.
  //  - on the rollers up to the lift-off point: the stinger arc, seated into the roller springs by the penetration
  //    that carries the contact pressure Te/R + w·cos φ;
  //  - free span: slope of the natural catenary at the current elevation, tan θ = sqrt(z(z + 2a))/a. When the pipe
  //    bears on the stinger tip the catenary leaves steeper than the tip tangent; the slope then starts at the tip
  //    tangent and approaches the catenary over the bending boundary layer λ = sqrt(EI/Te) (tensioned beam over a
  //    knife edge) instead of forming a kink;
  //  - on the seabed: horizontal, seated by the penetration w/kn that carries the weight.
  // Seating matters: nodes placed exactly on a contact surface are in or out of contact by rounding alone.
  const nn = n1 + n2 + n3 + 1, x = new Float64Array(nn), z = new Float64Array(nn);
  const kSeabed = o.seabed?.kn ?? 1e6;
  const catSlope = (zz) => Math.atan(Math.sqrt(zz * (zz + 2 * a)) / a);
  const kink = phiL < thTip ? 0 : phiL - catSlope(zL), lamB = Math.sqrt(pipe.EI / (H + w * zL));
  const arc = (p) => {
    const phi = p / R;
    const Rs = R - Math.min(1, phi / 0.02) * ((H + w * (zH - R * (1 - Math.cos(phi)))) / R + w * Math.cos(phi)) / rollerStiffness;
    return [Rs * Math.sin(phi), zH - R + Rs * Math.cos(phi)];
  };
  // The seabed end is pinned where the unstretched stinger arc + catenary + grounded length would end; the
  // tensioner end is free to move horizontally and takes up the stretch, so the march starts at the arc
  // coordinate p0 (on the flat tensioner bed when negative) that brings the last node to the pinned position.
  const n = nn - 1, xEnd = R * Math.sin(phiL) + a * Math.asinh(sig / a) + Lg;
  const march = (p0) => {
    let k = 0, p = p0, cx = Math.min(p0, 0), cz = zH, sFree = 0, sGround = -1;
    if (p0 > 0) [cx, cz] = arc(p0);
    x[0] = cx; z[0] = cz;
    for (const sec of sections) {
      const h = sec.length / sec.nElem;
      for (let i = 0; i < sec.nElem; i++) {
        let rem = h * (1 + (H + w * Math.max(cz, 0)) / pipe.EA);
        if (p < 0) { const d = Math.min(rem, -p); p += d; rem -= d; cx = p; }
        // on the arc the element is a chord, shorter than the arc it spans by the factor 1 − c²/(24R²)
        if (rem > 0 && p < sLift) { const gain = 1 + rem * rem / (24 * R * R), d = Math.min(rem, (sLift - p) / gain); p = Math.min(p + d * gain, sLift); rem -= d; [cx, cz] = arc(p); }
        const ds = rem / 8;
        for (let j = 0; j < 8 && rem > 0; j++) {
          if (sGround < 0) {
            const th0 = catSlope(cz) + kink * Math.exp(-sFree / lamB), zm = Math.max(cz - 0.5 * ds * Math.sin(th0), 0);
            const th = catSlope(zm) + kink * Math.exp(-(sFree + 0.5 * ds) / lamB);
            cx += ds * Math.cos(th); cz -= ds * Math.sin(th); sFree += ds;
            if (cz <= 0 || th <= 0) { cz = 0; sGround = 0; }
          } else { cx += ds; sGround += ds; cz = -Math.min(1, sGround / (0.02 * a)) * w / kSeabed; }
        }
        k++; x[k] = cx; z[k] = cz;
      }
    }
    return x[n];
  };
  let p0 = 0;
  for (let pass = 0; pass < 4; pass++) p0 += xEnd - march(p0);
  march(p0);
  const model = {
    sections,
    endA: { z: zH, angle: 0, Fx: -T }, endB: { x: xEnd, z: 0 },
    seabed: { z: 0, kn: 1e6, ...(o.seabed ?? {}) },
    supports: [{ xc: 0, zc: zH - R, R, k: rollerStiffness, angMin: Math.PI / 2 - thTip, angMax: Math.PI }],
    initial: { x, z }, rampLoads: false, rampSupports: false, nSteps: 1,
  };
  return { model, n, zL, thTip, hEl: lob / n1 };
}

/**
 * The rod.js model that sLay solves (sections, end conditions, seabed, stinger support and starting shape), for
 * callers that want to solve or vary it themselves with solveRodStatic.
 * @param {object} o sLay arguments
 * @returns {object} solveRodStatic model; `initial.x`, `initial.z` are the starting node coordinates [m]
 */
export function sLayRodModel(o) { return sLaySetup(o).model; }

/**
 * S-lay static configuration: pipe held by the tensioner at the top of a circular-arc stinger (roller bed of
 * radius R and arc length Ls), curving over the overbend, lifting off, passing through the sagbend and touching
 * down on the seabed.
 *
 * Method: rod.js beam model; the pipe end at the tensioner is clamped horizontal at the hinge elevation and
 * loaded by the tensioner force, the stinger is a frictionless rigid circular support active over its arc only
 * (the pipe may lift off before the last roller or bear on it), the seabed is a linear spring bed. The starting
 * shape is the stinger arc joined to a natural catenary with H = T − w·zHinge, pre-tensioned and seated on its
 * supports (see the comments in the code).
 * Reference: DNV-ST-F101 Sec. 10 (installation analysis); Bai & Bai, "Subsea Pipelines and Risers", Ch. 34;
 * analytic overbend strain ε = D/(2R) + T/EA.
 * Validity: static, in-plane, continuous roller support (discrete roller reactions not resolved), elastic pipe,
 * no stinger flexibility, no current.
 *
 * @param {object} o
 * @param {object} o.pipe section fields `{ EA, EI, w [N/m apparent weight, taken constant in air and water], … }`
 * @param {number} o.D steel outer diameter for bending strain [m]
 * @param {number} o.hingeHeight elevation of the stinger hinge / tensioner pipe axis above the seabed [m]
 * @param {number} o.stingerRadius R [m] @param {number} o.stingerLength arc length of the roller bed [m]
 * @param {number} o.tension tensioner force [N] (must exceed w·hingeHeight)
 * @param {number} [o.groundedLength] pipe modelled on the seabed [m] (default 35 % of the suspended length)
 * @param {object} [o.seabed] `{ kn [Pa] = 1e6, mu?, kt? }` @param {number} [o.rollerStiffness=1e8] [Pa]
 * @param {number} [o.nElem=260]
 * @returns {{ topTension: number, bottomTension: number, overbendCurvature: number, overbendStrain: number,
 *   overbendStrainAnalytic: number, sagbendCurvature: number, sagbendMoment: number, sagbendStrain: number,
 *   sSagbend: number, liftOffS: number, departureAngleDeg: number, bearsOnTip: boolean, stingerTipAngleDeg: number,
 *   tdpX: number, layBack: number, layBackFromHinge: number, suspendedLength: number, maxRollerReaction: number,
 *   rod: object }} tensions [N]; curvatures [1/m]; strains [-]; moment [N·m]; lengths/positions [m] (layBack =
 *   horizontal distance stinger tip → touchdown); angles [deg below horizontal]; maxRollerReaction [N/m]
 */
export function sLay(o) {
  const { model, n, zL, thTip, hEl } = sLaySetup(o);
  const { pipe, D, hingeHeight: zH, stingerRadius: R, tension: T } = o;
  const rod = solveRodStatic(model);
  let kOb = 0, kSag = 0, kLift = 0, kTdp = n, maxQ = 0;
  for (let i = 0; i <= n; i++) {
    if (rod.curvature[i] < rod.curvature[kOb]) kOb = i;
    if (rod.curvature[i] > rod.curvature[kSag]) kSag = i;
    if (rod.reaction[i] > 0 && rod.z[i] > 0.5 * zL) { kLift = i; maxQ = Math.max(maxQ, rod.reaction[i]); }
  }
  for (let i = kLift + 1; i <= n; i++) if (rod.reaction[i] > 0) { kTdp = i; break; }
  const xTip = R * Math.sin(thTip);
  const polar = Math.atan2(rod.z[kLift] - (zH - R), rod.x[kLift]);
  return {
    topTension: rod.Te[0], bottomTension: rod.Te[n],
    overbendCurvature: -rod.curvature[kOb], overbendStrain: -rod.curvature[kOb] * D / 2 + rod.Te[kOb] / pipe.EA,
    overbendStrainAnalytic: D / (2 * R) + T / pipe.EA,
    sagbendCurvature: rod.curvature[kSag], sagbendMoment: rod.moment[kSag],
    sagbendStrain: rod.curvature[kSag] * D / 2 + rod.Te[kSag] / pipe.EA, sSagbend: rod.s[kSag],
    liftOffS: rod.s[kLift], departureAngleDeg: -rod.angle[kLift] / DEG,
    bearsOnTip: (polar - (Math.PI / 2 - thTip)) * R < 1.5 * hEl, stingerTipAngleDeg: thTip / DEG,
    tdpX: rod.x[kTdp], layBack: rod.x[kTdp] - xTip, layBackFromHinge: rod.x[kTdp], suspendedLength: rod.s[kTdp] - rod.s[kLift],
    maxRollerReaction: maxQ, rod,
  };
}

/**
 * Minimum S-lay tension that keeps the sagbend bending strain within a limit (root search on sLay; the sagbend
 * strain falls monotonically with tension).
 * @param {object} o sLay arguments (tension ignored) plus `{ sagbendStrainLimit [-], tensionMin [N], tensionMax [N],
 *   tol? = 1e-3 [relative] }`
 * @returns {{ tension: number, result: object }} [N] and the sLay result at that tension
 */
export function minimumLayTension(o) {
  const { sagbendStrainLimit, tensionMin, tensionMax, tol = 1e-3 } = o ?? {};
  num(sagbendStrainLimit, 'sagbendStrainLimit', { min: 0, strictMin: true });
  num(tensionMin, 'tensionMin', { min: 0, strictMin: true }); num(tensionMax, 'tensionMax', { min: tensionMin, strictMin: true });
  const f = (T) => sLay({ ...o, tension: T }).sagbendStrain - sagbendStrainLimit;
  const T = brent(f, tensionMin, tensionMax, { tol, maxIter: 60 });
  return { tension: T, result: sLay({ ...o, tension: T }) };
}

/**
 * J-lay static configuration: pipe leaving a tower at a given angle from the vertical, with the natural
 * catenary matched to the tower angle (no overbend), bending stiffness and seabed springs.
 * Method: catenary sizing + rod model (riser.js solveScr with the top clamped at the tower angle).
 * Reference: DNV-ST-F101 Sec. 10; Bai & Bai Ch. 34. Validity: as solveScr.
 * @param {object} o `{ pipe, D [m], depth [m, tower exit height above the seabed], towerAngleDeg [deg from vertical],
 *   seabed?, nElem? = 240 }`
 * @returns {{ topTension: number, bottomTension: number, sagbendCurvature: number, sagbendStrain: number,
 *   sagbendMoment: number, layBack: number, suspendedLength: number, towerMoment: number, rod: object }}
 *   [N, N, 1/m, -, N·m, m, m, N·m]
 */
export function jLay({ pipe, D, depth, towerAngleDeg, seabed, nElem = 240 }) {
  num(D, 'D', { min: 0, strictMin: true }); num(towerAngleDeg, 'towerAngleDeg', { min: 0, strictMin: true, max: 80 });
  const r = solveScr({ pipe, h: depth, hangoffAngleDeg: towerAngleDeg, seabed, nElem, topCondition: { type: 'clamped', angleDeg: towerAngleDeg } });
  let kS = 0;
  for (let i = 0; i < r.s.length; i++) if (Math.abs(r.curvature[i]) > Math.abs(r.curvature[kS])) kS = i;
  return {
    topTension: r.topTension, bottomTension: r.bottomTension, sagbendCurvature: Math.abs(r.curvature[kS]),
    sagbendStrain: Math.abs(r.curvature[kS]) * D / 2 + r.Te[kS] / pipe.EA, sagbendMoment: Math.abs(r.moment[kS]),
    layBack: r.tdp ? r.tdp.x : NaN, suspendedLength: r.tdp ? r.tdp.s : NaN, towerMoment: Math.abs(r.moment[0]), rod: r,
  };
}

/**
 * Boundary-layer stiffened catenary (Plunkett 1967; Dixon & Rutledge 1968, "Stiffened catenary calculations in
 * pipeline laying problem", J. Eng. Ind. 90): the natural-catenary curvature κcat(σ) = a/(a² + σ²) (σ = arc length
 * from touchdown, a = H/w) multiplied by the end boundary layers (1 − e^(−σ/λ0)) at the touchdown, λ0 = sqrt(EI/H),
 * and (1 − e^(−(S−σ)/λt)) at a moment-free top, λt = sqrt(EI/Ttop). First-order composite expansion, valid for
 * λ ≪ S; used as an independent cross-check of the rod model.
 * @param {object} o `{ EI [N·m²], H [N], w [N/m], h [m top height above seabed], D? [m for strain], EA? [N],
 *   topFree? = true (moment-free top), nPoints? = 400 }`
 * @returns {{ lambda0: number, lambdaTop: number, suspendedLength: number, topTension: number,
 *   sagbendCurvature: number, sSagbend: number, sagbendMoment: number, sagbendStrain: number, sigma: Float64Array,
 *   curvature: Float64Array }} lengths [m]; tension [N]; curvature [1/m]; sSagbend = arc length from touchdown [m];
 *   strain [-] (NaN without D)
 */
export function stiffenedCatenary({ EI, H, w, h, D, EA, topFree = true, nPoints = 400 }) {
  num(EI, 'EI', { min: 0, strictMin: true }); num(H, 'H', { min: 0, strictMin: true });
  num(w, 'w', { min: 0, strictMin: true }); num(h, 'h', { min: 0, strictMin: true });
  const a = H / w, S = Math.sqrt(h * (h + 2 * a)), Tt = H + w * h;
  const l0 = Math.sqrt(EI / H), lt = Math.sqrt(EI / Tt);
  const sigma = new Float64Array(nPoints), curvature = new Float64Array(nPoints);
  let kM = 0;
  for (let i = 0; i < nPoints; i++) {
    const q = S * i / (nPoints - 1);
    sigma[i] = q;
    curvature[i] = (a / (a * a + q * q)) * (1 - Math.exp(-q / l0)) * (topFree ? 1 - Math.exp(-(S - q) / lt) : 1);
    if (curvature[i] > curvature[kM]) kM = i;
  }
  const Tsag = H + w * (Math.sqrt(a * a + sigma[kM] ** 2) - a);
  return {
    lambda0: l0, lambdaTop: lt, suspendedLength: S, topTension: Tt, sagbendCurvature: curvature[kM], sSagbend: sigma[kM],
    sagbendMoment: EI * curvature[kM], sagbendStrain: D ? curvature[kM] * D / 2 + (EA ? Tsag / EA : 0) : NaN, sigma, curvature,
  };
}

/* ------------------------------------------------------------------------------------------------------- */
/* Reel-lay                                                                                                  */

/**
 * Reel-lay bending: nominal bending strain on the reel, aligner and straightener, moment and elastic spring-back
 * of an elastic–perfectly-plastic thin-walled pipe, and ovalisation estimates.
 *
 * Strain: ε = (D/2)/(Rs + D/2) for a pipe wrapped on a surface of radius Rs (geometry).
 * Moment–curvature of a thin elastic–perfectly-plastic tube (first principles, r = mean radius, κy = fy/(E·r),
 * My = π·r²·t·fy): M/My = (2/π)·[(κ/κy)(θy − sinθy·cosθy) + 2·cosθy], sinθy = κy/κ (κ ≥ κy); residual curvature after
 * elastic unloading κres = κ − M/(E·I).
 * Ovalisation: elastic Brazier flattening ΔD/D = (1 − ν²)·r⁴·κ²/t² per side, f0 = 2·ΔD/D (Brazier 1927) — a lower
 * bound once the pipe yields; and the DNV-ST-F101 Sec. 5 bending ovalisation expression
 * f0' = [f0 + 0.030·(1 + D/(120t))·(2·ε·D/t)²]/(1 − pe/pc), evaluated here with the applied bending strain.
 * Validity: no strain hardening, no reel back-tension effect on ovality, single bend (cyclic accumulation over
 * reel/aligner/straightener not modelled).
 *
 * @param {object} o `{ D [m], t [m], reelRadius [m hub or current wrap surface radius], alignerRadius? [m],
 *   straightenerRadius? [m], fy [Pa], E? = 207e9, nu? = 0.3, f0? = 0.005 [-] initial ovality, pe? = 0, pc? [Pa] }`
 * @returns {{ reelStrain: number, alignerStrain: number, straightenerStrain: number, yieldCurvature: number,
 *   reelCurvature: number, reelMoment: number, residualCurvature: number, residualRadius: number,
 *   ovalityBrazier: number, ovalityF101: number }} strains [-]; curvatures [1/m]; moment [N·m]; radius [m];
 *   ovalities (Dmax − Dmin)/D [-]
 */
export function reelLay({ D, t, reelRadius, alignerRadius, straightenerRadius, fy, E = 207e9, nu = 0.3, f0 = 0.005, pe = 0, pc }) {
  num(D, 'D', { min: 0, strictMin: true }); num(t, 't', { min: 0, strictMin: true, max: D / 2 });
  num(reelRadius, 'reelRadius', { min: 0, strictMin: true }); num(fy, 'fy', { min: 0, strictMin: true });
  const strain = (Rs) => (Rs === undefined ? NaN : (D / 2) / (num(Rs, 'radius', { min: 0, strictMin: true }) + D / 2));
  const r = (D - t) / 2, I = Math.PI * (D ** 4 - (D - 2 * t) ** 4) / 64;
  const kap = 1 / (reelRadius + D / 2), ky = fy / (E * r), My = Math.PI * r * r * t * fy;
  let M;
  if (kap <= ky) M = My * kap / ky;
  else { const th = Math.asin(ky / kap); M = My * (2 / Math.PI) * ((kap / ky) * (th - Math.sin(th) * Math.cos(th)) + 2 * Math.cos(th)); }
  const kres = kap - M / (E * I);
  const eps = strain(reelRadius);
  const press = pc !== undefined ? 1 - pe / num(pc, 'pc', { min: 0, strictMin: true }) : 1;
  return {
    reelStrain: eps, alignerStrain: strain(alignerRadius), straightenerStrain: strain(straightenerRadius),
    yieldCurvature: ky, reelCurvature: kap, reelMoment: M, residualCurvature: kres, residualRadius: kres > 0 ? 1 / kres : Infinity,
    ovalityBrazier: f0 + 2 * (1 - nu * nu) * r ** 4 * kap * kap / (t * t),
    ovalityF101: (f0 + 0.030 * (1 + D / (120 * t)) * (2 * eps * D / t) ** 2) / press,
  };
}

/* ------------------------------------------------------------------------------------------------------- */
/* Cable / umbilical lay                                                                                     */

/**
 * Cable or umbilical lay over a chute or lay wheel: catenary tensions and geometry from one of bottom tension,
 * top tension or departure angle, chute bend-radius and sidewall-pressure checks, and slack.
 *
 * Catenary: Ttop = H + w·h, departure angle α from horizontal with cos α = H/Ttop, lay-back = a·asinh(S/a),
 * S = sqrt(h(h + 2a)), a = H/w; touchdown bend radius = a. Sidewall pressure on the chute SWP = T/R [N/m].
 * Slack (over-length) = payout speed / vessel speed − 1.
 * Reference: catenary statics (Irvine §1.2); DNV-RP-0360 §4 (cable installation: MBR, sidewall pressure);
 * Worzyk, "Submarine Power Cables", Ch. 7. Validity: no current, no bending stiffness, flat seabed.
 *
 * @param {object} o `{ w [N/m submerged weight], h [m chute departure height above seabed], H? | topTension? |
 *   departureAngleDeg? [deg from horizontal], chuteRadius? [m], mbr? [m], swpLimit? [N/m], payoutSpeed?, vesselSpeed?
 *   [m/s] }`
 * @returns {{ H: number, topTension: number, departureAngleDeg: number, layBack: number, suspendedLength: number,
 *   tdpBendRadius: number, sidewallPressure: number, ucChuteRadius: number, ucTdpRadius: number, ucSidewall: number,
 *   slack: number }} tensions [N]; angle [deg]; lengths [m]; SWP [N/m]; utilisations [-] (NaN when the limit or
 *   radius is not given); slack [-]
 */
export function cableLay({ w, h, H, topTension, departureAngleDeg, chuteRadius, mbr, swpLimit, payoutSpeed, vesselSpeed }) {
  num(w, 'w', { min: 0, strictMin: true }); num(h, 'h', { min: 0, strictMin: true });
  let Hh;
  if (H !== undefined) Hh = num(H, 'H', { min: 0, strictMin: true });
  else if (topTension !== undefined) Hh = num(topTension, 'topTension', { min: w * h, strictMin: true }) - w * h;
  else { const c = Math.cos(num(departureAngleDeg, 'departureAngleDeg', { min: 0, strictMin: true, max: 89.999 }) * DEG); Hh = w * h * c / (1 - c); }
  const cat = catenaryRiser({ h, w, H: Hh });
  const swp = chuteRadius !== undefined ? cat.topTension / num(chuteRadius, 'chuteRadius', { min: 0, strictMin: true }) : NaN;
  return {
    H: Hh, topTension: cat.topTension, departureAngleDeg: 90 - cat.hangoffAngleDeg, layBack: cat.horizontalSpan,
    suspendedLength: cat.suspendedLength, tdpBendRadius: cat.a, sidewallPressure: swp,
    ucChuteRadius: mbr !== undefined && chuteRadius !== undefined ? mbr / chuteRadius : NaN,
    ucTdpRadius: mbr !== undefined ? mbr / cat.a : NaN,
    ucSidewall: swpLimit !== undefined ? swp / swpLimit : NaN,
    slack: payoutSpeed !== undefined && vesselSpeed !== undefined ? payoutSpeed / num(vesselSpeed, 'vesselSpeed', { min: 0, strictMin: true }) - 1 : NaN,
  };
}

/**
 * Steady slack lay of a cable with zero bottom tension: the cable is a straight line at the critical incidence
 * angle α (from horizontal) where the normal component of weight balances normal drag,
 * w·cos α = ½ρ·Cd·d·(V·sin α)²  →  cos α = (−w + sqrt(w² + 4k²))/(2k), k = ½ρ·Cd·d·V²; top tension ≈ w·h
 * (tangential drag neglected). Reference: Zajac (1957), "Dynamics and kinematics of the laying and recovery of
 * submarine cable", Bell Syst. Tech. J. 36. Validity: steady speed, uniform cable, no current.
 * @param {object} o `{ w [N/m], d [m], V [m/s], h [m], Cd? = 1.2 [-], rho? = 1025 }`
 * @returns {{ angleDeg: number, topTension: number, layBack: number, suspendedLength: number }} [deg, N, m, m]
 */
export function cableSlackLay({ w, d, V, h, Cd = 1.2, rho = 1025 }) {
  num(w, 'w', { min: 0, strictMin: true }); num(d, 'd', { min: 0, strictMin: true }); num(V, 'V', { min: 0, strictMin: true });
  num(h, 'h', { min: 0, strictMin: true });
  const k = 0.5 * rho * Cd * d * V * V;
  const c = (-w + Math.sqrt(w * w + 4 * k * k)) / (2 * k), al = Math.acos(c);
  return { angleDeg: al / DEG, topTension: w * h, layBack: h / Math.tan(al), suspendedLength: h / Math.sin(al) };
}

/**
 * Abandonment and recovery: tension at the vessel in the A&R wire as the pipe (or cable) head is lowered or
 * raised, for a pipe catenary held at constant horizontal tension H. Head tension = H + wPipe·z (z = head height
 * above the seabed); vessel tension = head tension + wWire·(depth − z) (tension in a frictionless hanging line
 * varies by weight × elevation). Reference: catenary statics (Irvine §1.2); DNV-ST-F101 Sec. 10 (A&R).
 * Validity: static, no current, wire and pipe in one vertical plane, head fitting weight added as `headWeight`.
 * @param {object} o `{ H [N], wPipe [N/m], wWire [N/m], depth [m vessel sheave height above seabed],
 *   headHeights: number[] [m above seabed], headWeight? = 0 [N submerged] }`
 * @returns {Array<{ headHeight: number, headTension: number, vesselTension: number, wireAngleDeg: number,
 *   pipeSuspendedLength: number }>} [m, N, N, deg from horizontal at the head, m]
 */
export function abandonmentRecovery({ H, wPipe, wWire, depth, headHeights, headWeight = 0 }) {
  num(H, 'H', { min: 0, strictMin: true }); num(wPipe, 'wPipe', { min: 0, strictMin: true }); num(wWire, 'wWire', { min: 0 });
  num(depth, 'depth', { min: 0, strictMin: true }); num(headWeight, 'headWeight', { min: 0 });
  if (!Array.isArray(headHeights)) throw new RangeError('pipelay: argument headHeights must be an array');
  return headHeights.map((z, i) => {
    num(z, `headHeights[${i}]`, { min: 0, max: depth });
    const a = H / wPipe, S = Math.sqrt(z * (z + 2 * a));
    const V = wPipe * S + headWeight;                 // vertical force carried by the wire at the head
    const Th = Math.hypot(H, V);
    return { headHeight: z, headTension: Th, vesselTension: Th + wWire * (depth - z), wireAngleDeg: Math.atan2(V, H) / DEG, pipeSuspendedLength: S };
  });
}

/**
 * Capstan (Euler–Eytelwein) relation for a line bearing on a curved surface: T2 = T1·e^(μθ).
 * @param {object} o `{ T1 [N], mu [-], theta [rad wrap angle] }` @returns {number} T2 [N]
 */
export function capstan({ T1, mu, theta }) {
  num(T1, 'T1', { min: 0 }); num(mu, 'mu', { min: 0 }); num(theta, 'theta', { min: 0 });
  return T1 * Math.exp(mu * theta);
}

/**
 * Shore pull-in (or any pull along a route) tension build-up from the free tail to the winch:
 * straight segment T_out = T_in + w·L·(μ·cos β + sin β) (β = upward slope in the pull direction), bend
 * T_out = T_in·e^(μθ) (capstan). Reference: Coulomb friction and the capstan equation; DNV-RP-0360 §4 (pull-in).
 * Validity: quasi-static, no line bending stiffness or dynamic start-up factor (apply as `startFactor`).
 * @param {object} o `{ tailTension? = 0 [N], startFactor? = 1 [-], segments: Array<{ type: 'straight', length [m],
 *   w [N/m weight in the local medium], mu [-], slopeDeg? = 0 } | { type: 'bend', mu [-], angleDeg [deg] }> }`
 *   listed from the tail towards the winch
 * @returns {{ winchTension: number, tensions: number[] }} [N]; tension after each segment [N]
 */
export function shorePullIn({ tailTension = 0, startFactor = 1, segments }) {
  num(tailTension, 'tailTension', { min: 0 }); num(startFactor, 'startFactor', { min: 1 });
  if (!Array.isArray(segments) || segments.length === 0) throw new RangeError('pipelay: argument segments must be a non-empty array');
  let T = tailTension;
  const tensions = segments.map((sg, i) => {
    if (sg.type === 'bend') T = capstan({ T1: T, mu: num(sg.mu, `segments[${i}].mu`, { min: 0 }), theta: num(sg.angleDeg, `segments[${i}].angleDeg`, { min: 0 }) * DEG });
    else if (sg.type === 'straight') {
      const b = (sg.slopeDeg ?? 0) * DEG;
      T += num(sg.w, `segments[${i}].w`, { min: 0 }) * num(sg.length, `segments[${i}].length`, { min: 0 }) * (num(sg.mu, `segments[${i}].mu`, { min: 0 }) * Math.cos(b) + Math.sin(b));
      T = Math.max(T, 0);
    } else throw new RangeError(`pipelay: argument segments[${i}].type must be 'straight' or 'bend'`);
    return T;
  });
  return { winchTension: T * startFactor, tensions };
}

/**
 * Plough / trencher tow-force hook in the three-term form of Cathie & Wintgens (2001, OTC 13145):
 * F = Cw·W' + Cs·γ'·d³ + Cd·v·d² (interface friction + static soil resistance + rate effect). The coefficients are
 * site- and plough-specific and have no defaults here.
 * @param {object} o `{ Cw [-], Wsub [N plough submerged weight], Cs [-], gammaSub [N/m³], depth [m trench depth],
 *   Cd [N·s/m³ → force per (m/s·m²)], speed [m/s] }`
 * @returns {{ force: number, friction: number, staticSoil: number, rate: number }} [N]
 */
export function ploughTowForce({ Cw, Wsub, Cs, gammaSub, depth, Cd, speed }) {
  for (const [v, n] of [[Cw, 'Cw'], [Wsub, 'Wsub'], [Cs, 'Cs'], [gammaSub, 'gammaSub'], [depth, 'depth'], [Cd, 'Cd'], [speed, 'speed']]) num(v, n, { min: 0 });
  const friction = Cw * Wsub, staticSoil = Cs * gammaSub * depth ** 3, rate = Cd * speed * depth * depth;
  return { force: friction + staticSoil + rate, friction, staticSoil, rate };
}

/* ------------------------------------------------------------------------------------------------------- */
/* Quasi-dynamic tension and criteria                                                                        */

/**
 * Quasi-dynamic tension bounds for harmonic motion of the lay point (stinger tip, chute or tower top).
 *
 * Documented simple method (first principles, bounds rather than a prediction):
 *  - lower bound, line follows the motion quasi-statically: ΔT = w·a_v (catenary T = H + w·h) plus the inertia of the
 *    suspended line moved axially, m·Ls·ω²·a_t;
 *  - upper bound, "drag-locked": transverse drag prevents the catenary from changing shape within a cycle, so the
 *    axial component of the top motion is taken by elastic stretch, ΔT = (EA/Ls)·a_t.
 * a_v is the vertical amplitude, a_t the component of the motion amplitude along the line at the top
 * (a_t = |a_v·sin α + a_h·cos α|, α = top angle from horizontal). Pitch is included by the caller as
 * a_v = heave + arm·pitch. Validity: single harmonic, small motion; use rodHarmonicResponse for a resolved estimate.
 *
 * @param {object} o `{ staticTension [N], w [N/m], EA [N], m [kg/m incl. added mass in the axial direction],
 *   suspendedLength [m], topAngleDeg [deg from horizontal], heaveAmp [m], surgeAmp? = 0 [m], period [s] }`
 * @returns {{ axialAmp: number, lower: number, upper: number, dafLower: number, dafUpper: number,
 *   minTensionUpper: number }} amplitude [m]; dynamic tension amplitudes [N]; DAF = (Ts + ΔT)/Ts [-]; minimum
 *   tension with the upper-bound amplitude [N] (negative → compression/slack risk)
 */
export function layDynamicTension({ staticTension, w, EA, m, suspendedLength, topAngleDeg, heaveAmp, surgeAmp = 0, period }) {
  num(staticTension, 'staticTension', { min: 0, strictMin: true }); num(w, 'w', { min: 0 }); num(EA, 'EA', { min: 0, strictMin: true });
  num(m, 'm', { min: 0 }); num(suspendedLength, 'suspendedLength', { min: 0, strictMin: true });
  num(heaveAmp, 'heaveAmp', { min: 0 }); num(surgeAmp, 'surgeAmp', { min: 0 }); num(period, 'period', { min: 0, strictMin: true });
  const al = num(topAngleDeg, 'topAngleDeg', { min: 0, max: 90 }) * DEG, om = 2 * Math.PI / period;
  const at = Math.abs(heaveAmp * Math.sin(al) + surgeAmp * Math.cos(al));
  const lower = w * heaveAmp + m * suspendedLength * om * om * at;
  const upper = Math.max(lower, EA / suspendedLength * at);
  return { axialAmp: at, lower, upper, dafLower: 1 + lower / staticTension, dafUpper: 1 + upper / staticTension, minTensionUpper: staticTension - upper };
}

/**
 * Allowable sea-state search hook: the limiting dynamic amplification factor from tension limits,
 * DAF_lim = min(Tmax/Ts, 2 − Tmin/Ts) (maximum tension not above Tmax and minimum tension Ts·(2 − DAF) not below
 * Tmin), and — when a caller function DAF(Hs) is supplied — the largest Hs on the supplied grid that satisfies it.
 * @param {object} o `{ staticTension [N], maxTension [N], minTension? = 0 [N], hsGrid?: number[] [m],
 *   dafOfHs?: (Hs) => DAF }`
 * @returns {{ limitingDaf: number, allowableHs: number }} [-], [m] (NaN without the hook; 0 when no grid value passes)
 */
export function limitingSeaState({ staticTension, maxTension, minTension = 0, hsGrid, dafOfHs }) {
  num(staticTension, 'staticTension', { min: 0, strictMin: true }); num(maxTension, 'maxTension', { min: staticTension });
  num(minTension, 'minTension', { max: staticTension });
  const lim = Math.min(maxTension / staticTension, 2 - minTension / staticTension);
  let hs = NaN;
  if (Array.isArray(hsGrid) && typeof dafOfHs === 'function') { hs = 0; for (const v of hsGrid) if (dafOfHs(v) <= lim && v > hs) hs = v; }
  return { limitingDaf: lim, allowableHs: hs };
}

const OVERBEND_LIMITS = { X70: [0.00270, 0.00325], X65: [0.00250, 0.00305], X60: [0.00230, 0.00290], X52: [0.00205, 0.00260] };

/**
 * DNV-ST-F101 simplified laying criteria (Sec. 10 / former OS-F101 Sec. 13 H300).
 *
 * Overbend: total strain ≤ Criterion I (static loading) or Criterion II (static + dynamic) limit by grade —
 * X70 0.270 % / 0.325 %, X65 0.250 % / 0.305 %, X60 0.230 % / 0.290 %, X52 0.205 % / 0.260 % (documented defaults;
 * pass `strainLimit` to override). Sagbend and stinger tip: equivalent stress ≤ 0.87·fy with all load-effect factors
 * 1.0, σeq = sqrt(σl² + σh² − σl·σh), σl = N/A ± M·D/(2I) (N = Te − pe·Ae, empty pipe), σh = −pe·(D − t)/(2t).
 * @param {object} o `{ grade? 'X52' | 'X60' | 'X65' | 'X70', criterion? = 1 | 2, strainLimit? [-], overbendStrain [-],
 *   D [m], t [m], fy [Pa], sagbendTe [N effective tension], sagbendMoment [N·m], pe? = 0 [Pa], stressFactor? = 0.87 }`
 * @returns {{ overbendLimit: number, ucOverbend: number, sagbendStress: number, sagbendAllowable: number,
 *   ucSagbend: number }} [-], [-], [Pa], [Pa], [-]
 */
export function layCriteriaF101({ grade, criterion = 1, strainLimit, overbendStrain, D, t, fy, sagbendTe = 0, sagbendMoment = 0, pe = 0, stressFactor = 0.87 }) {
  num(D, 'D', { min: 0, strictMin: true }); num(t, 't', { min: 0, strictMin: true, max: D / 2 }); num(fy, 'fy', { min: 0, strictMin: true });
  num(overbendStrain, 'overbendStrain', { min: 0 }); num(pe, 'pe', { min: 0 });
  let lim = strainLimit;
  if (lim === undefined) {
    if (!OVERBEND_LIMITS[grade]) throw new RangeError("pipelay: argument grade must be 'X52', 'X60', 'X65' or 'X70' (or give strainLimit)");
    if (criterion !== 1 && criterion !== 2) throw new RangeError('pipelay: argument criterion must be 1 or 2');
    lim = OVERBEND_LIMITS[grade][criterion - 1];
  }
  num(lim, 'strainLimit', { min: 0, strictMin: true });
  const Di = D - 2 * t, A = Math.PI * (D * D - Di * Di) / 4, I = Math.PI * (D ** 4 - Di ** 4) / 64, Ae = Math.PI * D * D / 4;
  const N = sagbendTe - pe * Ae, sb = Math.abs(sagbendMoment) * D / (2 * I), sh = -pe * (D - t) / (2 * t);
  const eq = (sl) => Math.sqrt(sl * sl + sh * sh - sl * sh);
  const seq = Math.max(eq(N / A + sb), eq(N / A - sb));
  return { overbendLimit: lim, ucOverbend: overbendStrain / lim, sagbendStress: seq, sagbendAllowable: stressFactor * fy, ucSagbend: seq / (stressFactor * fy) };
}
