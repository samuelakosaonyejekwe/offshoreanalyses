/**
 * Slope stability: infinite slope, circular slip surfaces (ordinary method of slices and Bishop's
 * simplified method) with critical-circle search, wave-loaded submarine slopes, trench side slopes and
 * Newmark sliding-block displacement.
 *
 * SI units: m, N/m^3, Pa, rad, m/s^2, s. Two-dimensional sections; x horizontal, y up.
 *
 * References
 *   Fellenius W. (1936) Calculation of the stability of earth dams, 2nd Congress on Large Dams.
 *   Bishop A.W. (1955) The use of the slip circle in the stability analysis of slopes, Geotechnique 5(1).
 *   Taylor D.W. (1937) Stability of earth slopes, J. Boston Soc. Civ. Eng. 24.
 *   Skempton A.W. & DeLory F.A. (1957) Stability of natural slopes in London Clay, 4th ICSMFE.
 *   Newmark N.M. (1965) Effects of earthquakes on dams and embankments, Geotechnique 15(2).
 *   Henkel D.J. (1970) The role of waves in causing submarine landslides, Geotechnique 20(1).
 *   Duncan J.M., Wright S.G. & Brandon T.L. (2014) Soil Strength and Slope Stability, 2nd ed.
 */

const G0 = 9.80665;

function need(name, v, min = -Infinity, max = Infinity) {
  if (typeof v !== 'number' || Number.isNaN(v) || v < min || v > max) {
    throw new RangeError(`${name} must be a number in [${min}, ${max}] (got ${v})`);
  }
  return v;
}

/**
 * Infinite-slope factor of safety on a plane at depth z parallel to the surface (Skempton & DeLory 1957,
 * with pseudo-static seismic load).
 *   drained:    FS = [c′ + (γ z cos²β − kh γt z sin β cos β − u) tan φ′] / (γ z sin β cos β + kh γt z cos²β)
 *   undrained:  FS = su / (γ z sin β cos β + kh γt z cos²β)
 * γ is the unit weight that drives the slope (buoyant γ′ for a submerged slope without excess pore
 * pressure, bulk γ on land); γt the total unit weight for the inertia force (default γ). Pore pressure on
 * the plane: `u` directly, or `ru` (u = ru γ z), or seepage parallel to the slope with the phreatic
 * surface `hw` above the plane (u = γw hw cos²β).
 * @param {object} o
 * @param {number} o.beta slope angle [rad]
 * @param {number} o.z depth of the slip plane (vertical) [m]
 * @param {number} o.gamma driving unit weight [N/m^3]
 * @param {number} [o.gammaTotal=o.gamma] total unit weight for seismic inertia [N/m^3]
 * @param {'drained'|'undrained'} [o.mode='drained']
 * @param {number} [o.c=0] [Pa]  @param {number} [o.phi=0] [rad]  @param {number} [o.su] [Pa]
 * @param {number} [o.u] pore pressure on the plane [Pa]
 * @param {number} [o.ru] pore-pressure ratio [-]
 * @param {number} [o.hw] height of the phreatic surface above the plane (parallel seepage) [m]
 * @param {number} [o.gammaW=9810] [N/m^3]
 * @param {number} [o.kh=0] horizontal seismic coefficient [-]
 * @returns {{fs:number,tau:number,sigmaN:number,u:number,strength:number}} stresses on the plane [Pa]
 */
export function infiniteSlope({ beta, z, gamma, gammaTotal, mode = 'drained', c = 0, phi = 0, su, u, ru, hw, gammaW = 9810, kh = 0 }) {
  need('beta', beta, 0, Math.PI / 2 - 1e-9);
  need('z', z, 1e-12);
  need('gamma', gamma, 1e-12);
  need('kh', kh, 0, 1);
  const gt = gammaTotal === undefined ? gamma : gammaTotal;
  const sb = Math.sin(beta);
  const cb = Math.cos(beta);
  const tau = gamma * z * sb * cb + kh * gt * z * cb * cb;
  let uu = 0;
  if (u !== undefined) uu = u;
  else if (ru !== undefined) uu = ru * gamma * z;
  else if (hw !== undefined) uu = gammaW * need('hw', hw, 0) * cb * cb;
  const sigmaN = gamma * z * cb * cb - kh * gt * z * sb * cb;
  let strength;
  if (mode === 'undrained') strength = need('su', su, 0);
  else if (mode === 'drained') strength = c + Math.max(sigmaN - uu, 0) * Math.tan(phi);
  else throw new RangeError(`mode must be 'drained' or 'undrained' (got ${mode})`);
  return { fs: tau > 0 ? strength / tau : Infinity, tau, sigmaN, u: uu, strength };
}

/**
 * Yield (critical) acceleration of a sliding block from its static factor of safety (Newmark 1965):
 *   ky = (FS − 1) g sin β   for a planar slide parallel to a slope at β.
 * @param {number} fs static factor of safety [-]
 * @param {number} beta slope / thrust angle [rad]
 * @returns {number} yield acceleration [m/s^2] (0 if FS ≤ 1)
 */
export function yieldAcceleration(fs, beta) {
  need('fs', fs, 0);
  need('beta', beta, 0, Math.PI / 2);
  return Math.max(fs - 1, 0) * G0 * Math.sin(beta);
}

/**
 * Newmark (1965) rigid sliding-block displacement: the block slips whenever the ground acceleration
 * exceeds the yield acceleration and continues until the relative velocity returns to zero. One-way
 * (downslope) sliding; linear acceleration within each time step.
 * Check: a rectangular pulse of amplitude A > ky and duration t0 gives d = ½ (A − ky) t0² A / ky.
 * @param {object} o
 * @param {ArrayLike<number>} o.acc ground acceleration history [m/s^2]
 * @param {number} o.dt time step [s]
 * @param {number} o.ky yield acceleration [m/s^2], > 0
 * @param {boolean} [o.bothPolarities=false] also evaluate the reversed record and return the larger value
 * @returns {{displacement:number,maxVelocity:number,slipTime:number}} [m, m/s, s]
 */
export function newmarkDisplacement({ acc, dt, ky, bothPolarities = false }) {
  need('dt', dt, 1e-12);
  need('ky', ky, 1e-12);
  if (!acc || acc.length < 2) throw new RangeError('acc must contain at least two samples');
  const run = (sign) => {
    let v = 0;
    let d = 0;
    let vmax = 0;
    let slip = 0;
    // sub-stepping keeps the stick/slip transitions accurate for coarse records
    const sub = 10;
    const h = dt / sub;
    for (let i = 0; i < acc.length - 1; i++) {
      const a0 = sign * acc[i];
      const a1 = sign * acc[i + 1];
      for (let j = 0; j < sub; j++) {
        const ra = a0 + ((a1 - a0) * j) / sub - ky;
        const rb = a0 + ((a1 - a0) * (j + 1)) / sub - ky;
        if (v <= 0 && ra <= 0 && rb <= 0) { v = 0; continue; }
        let vn = v + 0.5 * (ra + rb) * h;
        if (v <= 0 && ra < 0) {
          // slip starts inside the sub-step, where the relative acceleration crosses zero
          const f = -ra / (rb - ra);
          vn = 0.5 * rb * (1 - f) * h;
          d += (vn * (1 - f) * h) / 3;
          slip += (1 - f) * h;
        } else if (vn < 0) {
          const f = v / (v - vn);
          d += 0.5 * v * f * h;
          slip += f * h;
          vn = 0;
        } else {
          d += 0.5 * (v + vn) * h;
          slip += h;
        }
        v = vn;
        if (v > vmax) vmax = v;
      }
    }
    return { displacement: d, maxVelocity: vmax, slipTime: slip };
  };
  const r = run(1);
  if (!bothPolarities) return r;
  const q = run(-1);
  return q.displacement > r.displacement ? q : r;
}

function makeGround(surface) {
  if (!Array.isArray(surface) || surface.length < 2) throw new RangeError('surface must be a polyline of at least two [x, y] points');
  for (let i = 1; i < surface.length; i++) if (!(surface[i][0] > surface[i - 1][0])) throw new RangeError('surface x-coordinates must be strictly increasing');
  const n = surface.length;
  return (x) => {
    if (x <= surface[0][0]) return surface[0][1];
    if (x >= surface[n - 1][0]) return surface[n - 1][1];
    let i = 1;
    while (surface[i][0] < x) i++;
    const [x0, y0] = surface[i - 1];
    const [x1, y1] = surface[i];
    return y0 + ((y1 - y0) * (x - x0)) / (x1 - x0);
  };
}

function normaliseLayers(soil) {
  const layers = Array.isArray(soil) ? [...soil] : [{ yTop: Infinity, ...soil }];
  layers.sort((a, b) => b.yTop - a.yTop);
  for (const l of layers) {
    need('layer.gamma', l.gamma, 0);
    need('layer.c', l.c === undefined ? 0 : l.c, 0);
    need('layer.phi', l.phi === undefined ? 0 : l.phi, 0, Math.PI / 2 - 1e-9);
  }
  return layers;
}

/**
 * Factor of safety of a circular slip surface by the ordinary method of slices (Fellenius) and Bishop's
 * simplified method:
 *   OMS:     FS = Σ [c l + (W cos α − u l) tan φ] / (Σ W sin α + Σ kh W e / R)
 *   Bishop:  FS = Σ { [c b + (W − u b) tan φ] / mα } / (Σ W sin α + Σ kh W e / R),  mα = cos α (1 + tan α tan φ / FS)
 * with α the base inclination (positive where the base dips in the direction of sliding), l = b / cos α,
 * e the lever arm of the slice centre about the circle centre (vertical distance). mα is limited to ≥ 0.2
 * to avoid the well-known numerical problem at steep negative base angles.
 * Soil: one homogeneous material {gamma, c, phi} or layers [{yTop, gamma, c, phi}] (a layer extends from
 * yTop down to the next layer). Pore pressure from `ru` or a piezometric line. For a fully submerged slope
 * use buoyant unit weights and no pore pressure. Tension cracks, external water loads and reinforcement are
 * not modelled.
 * @param {object} o
 * @param {[number,number][]} o.surface ground polyline, x increasing [m]
 * @param {{xc:number,yc:number,R:number}} o.circle slip circle [m]
 * @param {object|object[]} o.soil material(s): gamma [N/m^3], c [Pa], phi [rad], yTop [m]
 * @param {number} [o.nSlices=50]
 * @param {number} [o.kh=0] horizontal seismic coefficient [-]
 * @param {number} [o.ru] pore-pressure ratio u/(γh) [-]
 * @param {[number,number][]} [o.piezo] piezometric line [m]
 * @param {number} [o.gammaW=9810] [N/m^3]
 * @returns {{fsBishop:number,fsOrdinary:number,xEntry:number,xExit:number,drivingMoment:number,
 *   maxDepth:number,direction:number}|null} null when the circle does not cut a valid sliding mass;
 *   drivingMoment [N·m/m]; direction +1 = sliding towards +x
 */
export function circleStability({ surface, circle, soil, nSlices = 50, kh = 0, ru, piezo, gammaW = 9810 }) {
  const ground = makeGround(surface);
  const layers = normaliseLayers(soil);
  const pz = piezo ? makeGround(piezo) : null;
  return evalCircle(ground, surface, layers, circle.xc, circle.yc, circle.R, nSlices, kh, ru, pz, gammaW);
}

function evalCircle(ground, surface, layers, xc, yc, R, nSlices, kh, ru, pz, gammaW) {
  if (!(R > 0)) return null;
  const xMin = surface[0][0];
  const xMax = surface[surface.length - 1][0];
  const xl = Math.max(xc - R, xMin);
  const xr = Math.min(xc + R, xMax);
  if (!(xr > xl)) return null;
  const arc = (x) => yc - Math.sqrt(Math.max(R * R - (x - xc) * (x - xc), 0));
  const diff = (x) => ground(x) - arc(x);
  const m = 240;
  const dx = (xr - xl) / m;
  // widest contiguous run where the arc lies below the ground
  let best = null;
  let start = -1;
  for (let i = 0; i <= m + 1; i++) {
    const inside = i <= m && diff(xl + i * dx) > 0;
    if (inside && start < 0) start = i;
    if (!inside && start >= 0) {
      if (!best || i - 1 - start > best[1] - best[0]) best = [start, i - 1];
      start = -1;
    }
  }
  if (!best || best[1] - best[0] < 2) return null;
  const refine = (a, b) => {
    // a outside (diff ≤ 0), b inside
    for (let k = 0; k < 40; k++) { const mid = 0.5 * (a + b); if (diff(mid) > 0) b = mid; else a = mid; }
    return 0.5 * (a + b);
  };
  // both ends must daylight on the ground: an arc cut off by the edge of the section, or one whose
  // vertical tangent lies below ground (an unsupported vertical scarp), is not a valid slip surface
  if (best[0] === 0 || best[1] === m) return null;
  const xa = refine(xl + (best[0] - 1) * dx, xl + best[0] * dx);
  const xb = refine(xl + (best[1] + 1) * dx, xl + best[1] * dx);
  const b = (xb - xa) / nSlices;
  const W = new Float64Array(nSlices);
  const sinA = new Float64Array(nSlices);
  const cosA = new Float64Array(nSlices);
  const cc = new Float64Array(nSlices);
  const tp = new Float64Array(nSlices);
  const uu = new Float64Array(nSlices);
  let D = 0;
  let Ms = 0;
  let maxDepth = 0;
  for (let i = 0; i < nSlices; i++) {
    const xm = xa + (i + 0.5) * b;
    const yb = arc(xm);
    const yt = ground(xm);
    const h = Math.max(yt - yb, 0);
    if (h > maxDepth) maxDepth = h;
    let w = 0;
    let base = layers[layers.length - 1];
    for (let j = 0; j < layers.length; j++) {
      const top = Math.min(layers[j].yTop, yt);
      const bot = Math.max(j + 1 < layers.length ? layers[j + 1].yTop : -Infinity, yb);
      if (top > bot) w += layers[j].gamma * (top - bot);
      if (yb <= layers[j].yTop && (j + 1 >= layers.length || yb > layers[j + 1].yTop)) base = layers[j];
    }
    W[i] = w * b;
    sinA[i] = (xm - xc) / R;
    cosA[i] = Math.sqrt(Math.max(1 - sinA[i] * sinA[i], 1e-12));
    cc[i] = base.c || 0;
    tp[i] = Math.tan(base.phi || 0);
    uu[i] = ru !== undefined ? ru * w : pz ? Math.max(pz(xm) - yb, 0) * gammaW : 0;
    D += W[i] * sinA[i];
    Ms += kh * W[i] * (yc - (yb + h / 2)) / R;
  }
  if (Math.abs(D) < 1e-9) return null;
  // the arc centre lies above the mass; sliding is towards −x when Σ W sin α > 0 with this sign convention
  const dir = D > 0 ? -1 : 1;
  const s = D > 0 ? 1 : -1;
  const driving = Math.abs(D) + Ms;
  if (!(driving > 0)) return null;
  let num = 0;
  for (let i = 0; i < nSlices; i++) {
    const l = b / cosA[i];
    num += cc[i] * l + Math.max(W[i] * cosA[i] - uu[i] * l, 0) * tp[i];
  }
  const fsOrdinary = num / driving;
  let F = Math.max(fsOrdinary, 0.05);
  for (let it = 0; it < 200; it++) {
    let sum = 0;
    for (let i = 0; i < nSlices; i++) {
      const ma = Math.max(cosA[i] + (s * sinA[i] * tp[i]) / F, 0.2);
      sum += (cc[i] * b + (W[i] - uu[i] * b) * tp[i]) / ma;
    }
    const Fn = sum / driving;
    if (Math.abs(Fn - F) < 1e-9 * Math.max(Fn, 1)) { F = Fn; break; }
    F = Fn;
    if (!(F > 0)) return null;
  }
  return { fsBishop: F, fsOrdinary, xEntry: xa, xExit: xb, drivingMoment: driving * R, maxDepth, direction: dir };
}

/**
 * Search for the critical slip circle: coarse grid over centre and radius followed by a pattern-search
 * refinement. Returns the minimum factor of safety found (a local search cannot guarantee the global
 * minimum on irregular, strongly layered sections – widen the grid when in doubt).
 * @param {object} o surface, soil, nSlices, kh, ru, piezo, gammaW as circleStability, plus
 * @param {'bishop'|'ordinary'} [o.method='bishop']
 * @param {[number,number]} [o.xcRange] centre search range in x [m] (default: extent of the surface)
 * @param {[number,number]} [o.ycRange] centre search range in y [m] (default: 0.05 H to 3 H above the crest)
 * @param {number} [o.yRigid] elevation of a rigid stratum the circle may not penetrate [m]
 *   (default: one slope height below the lowest ground point)
 * @param {number} [o.minDepth] minimum depth of the sliding mass [m] (default 2 % of the slope height)
 * @param {number} [o.grid=14] grid points per dimension
 * @returns {{fs:number,circle:{xc:number,yc:number,R:number},fsBishop:number,fsOrdinary:number,xEntry:number,
 *   xExit:number,evaluations:number}}
 */
export function criticalCircle(o) {
  const { surface, soil, nSlices = 40, kh = 0, ru, piezo, gammaW = 9810, method = 'bishop', grid = 14 } = o;
  const ground = makeGround(surface);
  const layers = normaliseLayers(soil);
  const pz = piezo ? makeGround(piezo) : null;
  let yLo = Infinity;
  let yHi = -Infinity;
  for (const p of surface) { yLo = Math.min(yLo, p[1]); yHi = Math.max(yHi, p[1]); }
  const H = yHi - yLo;
  if (!(H > 0)) throw new RangeError('surface must have relief (a slope)');
  const xr = o.xcRange || [surface[0][0], surface[surface.length - 1][0]];
  const yr = o.ycRange || [yHi + 0.05 * H, yHi + 3 * H];
  const yRigid = o.yRigid !== undefined ? o.yRigid : yLo - H;
  const minDepth = o.minDepth !== undefined ? o.minDepth : 0.02 * H;
  let evals = 0;
  const key = method === 'ordinary' ? 'fsOrdinary' : 'fsBishop';
  const f = (xc, yc, R) => {
    if (yc - R < yRigid - 1e-9) return null;
    evals++;
    const r = evalCircle(ground, surface, layers, xc, yc, R, nSlices, kh, ru, pz, gammaW);
    if (!r || r.maxDepth < minDepth) return null;
    return r;
  };
  let best = null;
  const consider = (xc, yc, R) => {
    const r = f(xc, yc, R);
    if (r && (!best || r[key] < best.r[key])) { best = { r, xc, yc, R }; return true; }
    return false;
  };
  for (let i = 0; i < grid; i++) {
    const xc = xr[0] + ((xr[1] - xr[0]) * i) / (grid - 1);
    for (let j = 0; j < grid; j++) {
      const yc = yr[0] + ((yr[1] - yr[0]) * j) / (grid - 1);
      const Rmin = Math.max(yc - yHi, 1e-6);
      const Rmax = yc - yRigid;
      for (let k = 1; k <= grid + 4; k++) consider(xc, yc, Rmin + ((Rmax - Rmin) * k) / (grid + 4));
    }
  }
  if (!best) throw new Error('no valid slip circle found in the search domain');
  let sx = (xr[1] - xr[0]) / (grid - 1);
  let sy = (yr[1] - yr[0]) / (grid - 1);
  let sr = 0.5 * H / grid + (yr[1] - yLo) / (2 * grid);
  for (let it = 0; it < 60 && Math.max(sx, sy, sr) > 1e-4 * H; it++) {
    let improved = false;
    for (const [ax, ay, ar] of [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1], [0, 1, 1], [0, -1, -1], [1, 0, 1], [-1, 0, -1]]) {
      if (consider(best.xc + ax * sx, best.yc + ay * sy, best.R + ar * sr)) improved = true;
    }
    if (!improved) { sx *= 0.5; sy *= 0.5; sr *= 0.5; }
  }
  return { fs: best.r[key], circle: { xc: best.xc, yc: best.yc, R: best.R }, fsBishop: best.r.fsBishop, fsOrdinary: best.r.fsOrdinary, xEntry: best.r.xEntry, xExit: best.r.xExit, evaluations: evals };
}

/**
 * Side-slope stability of a dredged trench or rubble-mound slope: builds the section (crest plateau,
 * slope at `slopeAngle`, floor) and runs the critical-circle search. For cohesionless material the
 * infinite-slope value tan φ / tan β (dry or fully submerged, no seepage) is returned alongside.
 * Use buoyant unit weight for a submerged trench.
 * @param {object} o
 * @param {number} o.depth trench depth / slope height [m]
 * @param {number} o.slopeAngle side-slope angle from horizontal [rad]
 * @param {object|object[]} o.soil {gamma, c, phi} or layers with yTop measured with the crest at y = 0
 * @param {number} [o.kh=0]  @param {number} [o.ru]
 * @param {number} [o.yRigid] [m]
 * @returns {{fs:number,circle:object,fsInfinite:number|undefined}}
 */
export function trenchSideSlope({ depth, slopeAngle, soil, kh = 0, ru, yRigid }) {
  need('depth', depth, 1e-9);
  need('slopeAngle', slopeAngle, 1e-6, Math.PI / 2 - 1e-6);
  const run = depth / Math.tan(slopeAngle);
  // sliding towards +x: crest on the left, floor on the right
  const surface = [[-3 * depth - run, 0], [0, 0], [run, -depth], [run + 3 * depth + run, -depth]];
  const r = criticalCircle({ surface, soil, kh, ru, yRigid, xcRange: [-0.5 * run - depth, 1.5 * run + depth] });
  let fsInfinite;
  if (!Array.isArray(soil) && !(soil.c > 0) && soil.phi > 0) {
    fsInfinite = infiniteSlope({ beta: slopeAngle, z: 1, gamma: soil.gamma, phi: soil.phi, ru, kh }).fs;
  }
  return { fs: r.fs, circle: r.circle, fsInfinite };
}

/**
 * Undrained stability of a gently sloping seabed under wave loading (Henkel 1970 concept, elastic stress
 * field): at depth z the downslope gravity shear γ′ z sin β cos β is augmented by the wave-induced shear
 * stress amplitude in a deep homogeneous seabed, τw = p0 k z e^(−kz), p0 = γw H / (2 cosh kd):
 *     FS(z) = su(z) / (γ′ z sin β cos β + τw(z)).
 * The minimum over depth is returned. Valid for small slopes and seabed thickness larger than about half
 * a wavelength; for layered or shallow beds use a seabed-response analysis.
 * @param {object} o
 * @param {number} o.beta seabed slope [rad]
 * @param {number} o.H wave height [m]  @param {number} o.k wave number [rad/m]  @param {number} o.d water depth [m]
 * @param {number} o.gammaSub γ′ [N/m^3]
 * @param {number} o.su0 [Pa]  @param {number} [o.suGrad=0] [Pa/m]
 * @param {number} [o.zMax] depth range examined [m] (default one wavelength)
 * @param {number} [o.gammaW=10055] [N/m^3]
 * @param {number} [o.n=200]
 * @returns {{fsMin:number,zCritical:number,p0:number,z:number[],fs:number[]}}
 */
export function submarineSlopeWave({ beta, H, k, d, gammaSub, su0, suGrad = 0, zMax, gammaW = 10055, n = 200 }) {
  need('beta', beta, 0, Math.PI / 4);
  need('H', H, 0);
  need('k', k, 1e-12);
  need('d', d, 1e-9);
  need('gammaSub', gammaSub, 0);
  need('su0', su0, 0);
  const p0 = (gammaW * H) / (2 * Math.cosh(k * d));
  const zm = zMax !== undefined ? zMax : (2 * Math.PI) / k;
  const z = [];
  const fs = [];
  let fsMin = Infinity;
  let zCritical = 0;
  for (let i = 1; i <= n; i++) {
    const zi = (zm * i) / n;
    const tau = gammaSub * zi * Math.sin(beta) * Math.cos(beta) + p0 * k * zi * Math.exp(-k * zi);
    const v = tau > 0 ? (su0 + suGrad * zi) / tau : Infinity;
    z.push(zi);
    fs.push(v);
    if (v < fsMin) { fsMin = v; zCritical = zi; }
  }
  return { fsMin, zCritical, p0, z, fs };
}
