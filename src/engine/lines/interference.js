/**
 * Line / riser interference under steady current: quasi-static current deflection of a hanging or
 * top-tensioned line, wake shielding of a downstream cylinder (Huse and Blevins wake models), clearance between
 * two neighbouring lines along their length, collision screening and an impact-velocity estimate.
 *
 * Deflection model: near-vertical tensioned string (no bending stiffness, small angles). With s measured upward
 * from the lower end, effective tension Te(s) = Tbottom + w s and in-line drag per unit length
 * f(s) = 0.5 rho Cd d U(s)|U(s)|, the lateral equilibrium is d/ds (Te dy/ds) = -f.
 *
 * References:
 *  - Huse, E. (1993) "Interaction in deep-sea riser arrays", OTC 7237 (wake velocity deficit); also quoted in
 *    DNV-RP-F203 Riser Interference.
 *  - Blevins, R.D. (2005) "Forces on and stability of a cylinder in a wake", J. Offshore Mech. Arct. Eng. 127.
 *  - Sparks, C.P. (2007) Fundamentals of Marine Riser Mechanics, ch. 6-7 (tensioned-string riser under current).
 */

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}

function speedAt(current, s) {
  if (typeof current === 'number') return current;
  const { s: xs, speed: us } = current;
  const n = xs.length;
  if (s <= xs[0]) return us[0];
  if (s >= xs[n - 1]) return us[n - 1];
  let k = 0;
  while (s > xs[k + 1]) k++;
  return us[k] + (us[k + 1] - us[k]) * (s - xs[k]) / (xs[k + 1] - xs[k]);
}

function checkCurrent(current) {
  if (typeof current === 'number') { need(Number.isFinite(current), 'current must be finite'); return; }
  need(current && Array.isArray(current.s) && Array.isArray(current.speed) && current.s.length === current.speed.length && current.s.length >= 1,
    'current must be a speed [m/s] or { s: [], speed: [] } with equal lengths');
  for (let i = 1; i < current.s.length; i++) need(current.s[i] > current.s[i - 1], 'current.s must be ascending');
}

/**
 * Steady-current deflection of a near-vertical line or riser (quasi-static, drag only).
 *
 * mode 'top-tensioned': pinned at both ends, Te(s) = topTension - w (L - s);
 *   y(s) = integral_0^s (c - F)/Te ds with F(s) = integral_0^s f ds and c fixed by y(L) = 0.
 * mode 'free-hanging': pinned at the top, free lower end carrying endWeight (submerged) and an optional lateral
 *   end load; Te(s) = endWeight + w s; Te dy/ds = -(F + endForce), y(L) = 0.
 * Closed-form checks: constant tension -> parabola with y_mid = f L^2 / (8 T); free-hanging uniform line in
 * uniform current -> straight at tan(angle) = f / w.
 * Reference: Sparks (2007) ch. 6-7. Validity: small angles (check maxAngle, < ~0.2 rad), Te > 0 everywhere, no
 * bending stiffness, no lift, drag evaluated on the undeflected (vertical) line with the normal velocity = U.
 *
 * @param {object} o
 * @param {number} o.L length [m]
 * @param {number} o.w submerged weight per unit length [N/m] (>= 0)
 * @param {number} o.d drag diameter [m] @param {number} o.Cd drag coefficient [-]
 * @param {number|{s:number[],speed:number[]}} o.current speed [m/s], uniform or tabulated against height s above
 *   the lower end [m]
 * @param {'top-tensioned'|'free-hanging'} [o.mode='top-tensioned']
 * @param {number} [o.topTension] effective tension at the top [N] (top-tensioned)
 * @param {number} [o.endWeight=0] submerged weight hanging at the lower end [N] (free-hanging)
 * @param {number} [o.endForce=0] lateral current force on the end body [N] (free-hanging)
 * @param {number} [o.rho=1025] [kg/m^3] @param {number} [o.n=200] number of intervals
 * @param {ArrayLike<number>} [o.dragFactor] multiplier on f at each of the n + 1 stations (wake shielding) [-]
 * @returns {{s:Float64Array,y:Float64Array,slope:Float64Array,Te:Float64Array,f:Float64Array,yMax:number,
 *   sAtMax:number,maxAngle:number,topAngle:number,bottomAngle:number,totalDrag:number}} stations [m],
 *   downstream deflection [m], dy/ds [-], effective tension [N], drag per length [N/m], largest deflection [m]
 *   and its location, largest |angle| from vertical [rad], end angles [rad], total drag [N]
 */
export function currentDeflection(o) {
  const { L, w, d, Cd } = o;
  const mode = o.mode ?? 'top-tensioned';
  const rho = o.rho ?? 1025, n = o.n ?? 200;
  need(Number.isFinite(L) && L > 0, 'L must be > 0');
  need(Number.isFinite(w) && w >= 0, 'w must be >= 0');
  need(Number.isFinite(d) && d > 0, 'd must be > 0');
  need(Number.isFinite(Cd) && Cd >= 0, 'Cd must be >= 0');
  need(Number.isInteger(n) && n >= 2, 'n must be an integer >= 2');
  need(mode === 'top-tensioned' || mode === 'free-hanging', "mode must be 'top-tensioned' or 'free-hanging'");
  checkCurrent(o.current);
  if (o.dragFactor) need(o.dragFactor.length === n + 1, 'dragFactor must have n + 1 entries');
  let Tb;
  if (mode === 'top-tensioned') {
    need(Number.isFinite(o.topTension), 'topTension must be given for a top-tensioned line');
    Tb = o.topTension - w * L;
    need(Tb > 0, 'topTension must exceed the submerged weight w L (positive tension at the lower end)');
  } else {
    Tb = o.endWeight ?? 0;
    need(Tb >= 0, 'endWeight must be >= 0');
    need(Tb > 0 || w > 0, 'a free-hanging line needs w > 0 or endWeight > 0');
  }
  const Fend = mode === 'free-hanging' ? (o.endForce ?? 0) : 0;
  const h = L / n;
  const s = new Float64Array(n + 1), y = new Float64Array(n + 1), slope = new Float64Array(n + 1);
  const Te = new Float64Array(n + 1), f = new Float64Array(n + 1), F = new Float64Array(n + 1);
  for (let i = 0; i <= n; i++) {
    s[i] = i * h;
    const U = speedAt(o.current, s[i]);
    f[i] = 0.5 * rho * Cd * d * U * Math.abs(U) * (o.dragFactor ? o.dragFactor[i] : 1);
    Te[i] = Tb + w * s[i];
    if (i > 0) F[i] = F[i - 1] + 0.5 * h * (f[i] + f[i - 1]);
  }
  // g = (F + Fend)/Te with the finite limit f(0)/w at a tension-free lower end
  const g = new Float64Array(n + 1), inv = new Float64Array(n + 1);
  for (let i = 0; i <= n; i++) {
    if (Te[i] > 0) { g[i] = (F[i] + Fend) / Te[i]; inv[i] = 1 / Te[i]; } else { g[i] = Fend === 0 ? f[0] / w : Infinity; inv[i] = 0; }
  }
  need(Number.isFinite(g[0]), 'endForce needs a non-zero tension at the lower end (endWeight > 0)');
  let c = 0;
  if (mode === 'top-tensioned') {
    let a = 0, b = 0;
    for (let i = 1; i <= n; i++) { a += 0.5 * h * (g[i] + g[i - 1]); b += 0.5 * h * (inv[i] + inv[i - 1]); }
    c = a / b;
  }
  for (let i = 0; i <= n; i++) slope[i] = c * inv[i] - g[i];
  if (mode === 'top-tensioned') {
    for (let i = 1; i <= n; i++) y[i] = y[i - 1] + 0.5 * h * (slope[i] + slope[i - 1]);
    y[n] = 0;
  } else {
    for (let i = n - 1; i >= 0; i--) y[i] = y[i + 1] - 0.5 * h * (slope[i] + slope[i + 1]);
  }
  let yMax = 0, sAt = 0, aMax = 0;
  for (let i = 0; i <= n; i++) {
    if (Math.abs(y[i]) > Math.abs(yMax)) { yMax = y[i]; sAt = s[i]; }
    aMax = Math.max(aMax, Math.abs(Math.atan(slope[i])));
  }
  return {
    s, y, slope, Te, f, yMax, sAtMax: sAt, maxAngle: aMax, topAngle: Math.atan(slope[n]), bottomAngle: Math.atan(slope[0]),
    totalDrag: F[n] + Fend,
  };
}

/**
 * Huse wake model: velocity deficit behind a cylinder in steady flow,
 *   u(x, y) = k2 U0 sqrt(Cd D / xs) exp(-0.693 (y / b)^2),  b = k1 sqrt(Cd D xs),  xs = x + 4 D / Cd,
 * with k1 = 0.25 and k2 = 1.0 for a smooth cylinder (Huse 1993, OTC 7237; DNV-RP-F203). xs places a virtual wake
 * origin upstream so that the model is usable close to the cylinder. The local in-line velocity is U0 - u and
 * the drag on a downstream cylinder referred to the free-stream velocity is Cd_down (1 - u/U0)^2.
 * Valid for x > 0 (downstream), far-wake form; k1, k2 may be tuned to data.
 * @param {object} o
 * @param {number} o.x in-line distance downstream of the upstream cylinder centre [m] (> 0)
 * @param {number} [o.y=0] transverse distance from the wake centreline [m]
 * @param {number} o.D upstream cylinder diameter [m] @param {number} o.Cd upstream drag coefficient [-]
 * @param {number} [o.U0=1] free-stream speed [m/s]
 * @param {number} [o.k1=0.25] @param {number} [o.k2=1.0] model constants [-]
 * @returns {{deficit:number,velocity:number,halfWidth:number,xs:number,dragFactor:number}} velocity deficit u
 *   [m/s], local velocity U0 - u [m/s], wake half-width b [m], xs [m], (1 - u/U0)^2 [-]
 */
export function huseWake({ x, y = 0, D, Cd, U0 = 1, k1 = 0.25, k2 = 1.0 }) {
  need(Number.isFinite(x) && x > 0, 'x must be > 0 (downstream of the cylinder)');
  need(Number.isFinite(y), 'y must be finite');
  need(Number.isFinite(D) && D > 0, 'D must be > 0');
  need(Number.isFinite(Cd) && Cd > 0, 'Cd must be > 0');
  need(k1 > 0 && k2 > 0, 'k1 and k2 must be > 0');
  const xs = x + 4 * D / Cd;
  const b = k1 * Math.sqrt(Cd * D * xs);
  let ratio = k2 * Math.sqrt(Cd * D / xs) * Math.exp(-0.693 * (y / b) * (y / b));
  if (ratio > 1) ratio = 1;
  return { deficit: ratio * U0, velocity: (1 - ratio) * U0, halfWidth: b, xs, dragFactor: (1 - ratio) * (1 - ratio) };
}

/**
 * Blevins wake model: drag coefficient of a cylinder in the wake of an upstream cylinder, referred to the
 * free-stream velocity,
 *   Cd(x, y) = Cd0 [1 - a1 sqrt(Cd_up D_up / x) exp(-a2 y^2 / (Cd_up D_up x))]^2,
 * with a1 = 1.0 and a2 = 4.5 (Blevins 2005, JOMAE 127). The wake lift force of that paper is not included.
 * Valid for x greater than a few diameters; a1, a2 may be tuned to data.
 * @param {object} o
 * @param {number} o.x in-line spacing between cylinder centres [m] (> 0) @param {number} [o.y=0] transverse offset [m]
 * @param {number} o.D upstream diameter [m] @param {number} o.Cd upstream drag coefficient [-]
 * @param {number} [o.Cd0=o.Cd] free-stream drag coefficient of the downstream cylinder [-]
 * @param {number} [o.a1=1.0] @param {number} [o.a2=4.5] model constants [-]
 * @returns {{Cd:number,dragFactor:number}} downstream drag coefficient [-] and Cd / Cd0 [-]
 */
export function blevinsWakeDrag({ x, y = 0, D, Cd, Cd0, a1 = 1.0, a2 = 4.5 }) {
  need(Number.isFinite(x) && x > 0, 'x must be > 0');
  need(Number.isFinite(y), 'y must be finite');
  need(Number.isFinite(D) && D > 0, 'D must be > 0');
  need(Number.isFinite(Cd) && Cd > 0, 'Cd must be > 0');
  need(a1 > 0 && a2 > 0, 'a1 and a2 must be > 0');
  let r = a1 * Math.sqrt(Cd * D / x) * Math.exp(-a2 * y * y / (Cd * D * x));
  if (r > 1) r = 1;
  const fac = (1 - r) * (1 - r);
  return { Cd: (Cd0 ?? Cd) * fac, dragFactor: fac };
}

/**
 * Clearance between two neighbouring near-vertical lines / risers in a steady current flowing from the
 * upstream line towards the downstream one. The upstream line deflects in the free stream; the drag on the
 * downstream line is reduced by wake shielding evaluated at the local centre-to-centre distance, which itself
 * depends on both deflections (fixed-point iteration). Both lines must have the same length and share the same
 * stations (same end elevations).
 * Reference: DNV-RP-F203 (riser interference: wake models and clearance assessment); Huse (1993); Blevins (2005).
 * Validity: as currentDeflection; in-line arrangement with a constant transverse offset; no wake-induced lift,
 * no VIV drag amplification (apply it through Cd).
 * @param {object} o
 * @param {object} o.upstream line definition for currentDeflection (without current / n / dragFactor)
 * @param {object} o.downstream line definition (same L)
 * @param {number|{s:number[],speed:number[]}} o.current current as in currentDeflection
 * @param {number|ArrayLike<number>} o.spacing undeflected centre-to-centre in-line spacing [m] (one value or
 *   n + 1 values along the length)
 * @param {number} [o.transverseOffset=0] transverse offset of the downstream line from the wake centreline [m]
 * @param {'huse'|'blevins'|'none'} [o.wakeModel='huse']
 * @param {object} [o.wake] model constants ({ k1, k2 } or { a1, a2 })
 * @param {number} [o.n=200] intervals @param {number} [o.rho=1025] [kg/m^3]
 * @param {number} [o.maxIter=50] @param {number} [o.tol=1e-6] convergence on deflection [m]
 * @returns {{s:Float64Array,clearance:Float64Array,distance:Float64Array,dragFactor:Float64Array,
 *   upstream:object,downstream:object,minClearance:number,sAtMin:number,collision:boolean,iterations:number,
 *   converged:boolean}} surface-to-surface clearance and centre distance along the length [m]; shielding factor
 *   on the downstream drag [-]; the two deflection results; collision = minClearance <= 0
 */
export function clearanceUnderCurrent(o) {
  const n = o.n ?? 200, rho = o.rho ?? 1025;
  need(o.upstream && o.downstream, 'upstream and downstream line definitions are required');
  need(o.upstream.L === o.downstream.L, 'upstream.L and downstream.L must be equal');
  const model = o.wakeModel ?? 'huse';
  need(['huse', 'blevins', 'none'].includes(model), "wakeModel must be 'huse', 'blevins' or 'none'");
  const sp = new Float64Array(n + 1);
  if (typeof o.spacing === 'number') sp.fill(o.spacing);
  else { need(o.spacing && o.spacing.length === n + 1, 'spacing must be a number or n + 1 values'); sp.set(o.spacing); }
  for (let i = 0; i <= n; i++) need(Number.isFinite(sp[i]) && sp[i] > 0, 'spacing must be > 0');
  const yT = o.transverseOffset ?? 0;
  const up = currentDeflection({ ...o.upstream, current: o.current, n, rho });
  const fac = new Float64Array(n + 1).fill(1);
  let dn = currentDeflection({ ...o.downstream, current: o.current, n, rho });
  const dist = new Float64Array(n + 1);
  const rMin = 0.5 * (o.upstream.d + o.downstream.d);
  let it = 0, converged = model === 'none';
  const maxIter = o.maxIter ?? 50, tol = o.tol ?? 1e-6;
  for (; it < maxIter && !converged; it++) {
    for (let i = 0; i <= n; i++) {
      // shielding is evaluated no closer than contact; beyond contact the lines are flagged as colliding
      const x = Math.max(sp[i] + dn.y[i] - up.y[i], rMin);
      const target = model === 'huse'
        ? huseWake({ x, y: yT, D: o.upstream.d, Cd: o.upstream.Cd, ...(o.wake ?? {}) }).dragFactor
        : blevinsWakeDrag({ x, y: yT, D: o.upstream.d, Cd: o.upstream.Cd, ...(o.wake ?? {}) }).dragFactor;
      fac[i] = 0.5 * (fac[i] + target); // under-relaxed
    }
    const nx = currentDeflection({ ...o.downstream, current: o.current, n, rho, dragFactor: fac });
    let ch = 0;
    for (let i = 0; i <= n; i++) ch = Math.max(ch, Math.abs(nx.y[i] - dn.y[i]));
    dn = nx;
    if (ch < tol) converged = true;
  }
  const clr = new Float64Array(n + 1);
  let mn = Infinity, sAt = 0;
  for (let i = 0; i <= n; i++) {
    dist[i] = Math.hypot(sp[i] + dn.y[i] - up.y[i], yT);
    clr[i] = dist[i] - rMin;
    if (clr[i] < mn) { mn = clr[i]; sAt = up.s[i]; }
  }
  return {
    s: up.s, clearance: clr, distance: dist, dragFactor: fac, upstream: up, downstream: dn,
    minClearance: mn, sAtMin: sAt, collision: mn <= 0, iterations: it, converged,
  };
}

/**
 * Impact-velocity estimate for collision screening, from the work-energy theorem: a line section of mass m per
 * unit length (including added mass), initially at rest, is driven across a gap g by a net lateral load per unit
 * length f0 that is opposed by a linear restoring stiffness k per unit length:
 *   0.5 m v^2 = f0 g - 0.5 k g^2   ->   v = sqrt((2 f0 g - k g^2) / m)   (v = 0 when the gap is not closed).
 * Fluid damping is neglected, so v is an upper bound for this idealisation; it is a first-principles screening
 * estimate, not a code formula. For a tensioned string of span L at mid-span k ~ 8 T / L^2 (from
 * y_mid = f L^2 / (8 T)).
 * @param {object} o
 * @param {number} o.force net lateral driving load per unit length f0 [N/m]
 * @param {number} o.gap distance travelled to contact g [m]
 * @param {number} o.mass mass plus added mass per unit length [kg/m]
 * @param {number} [o.stiffness=0] lateral restoring stiffness per unit length k [N/m^2]
 * @returns {{velocity:number,energy:number,reaches:boolean}} impact velocity [m/s], kinetic energy per unit
 *   length at impact [J/m], whether the gap is closed
 */
export function impactVelocityEstimate({ force, gap, mass, stiffness = 0 }) {
  need(Number.isFinite(force) && force >= 0, 'force must be >= 0');
  need(Number.isFinite(gap) && gap >= 0, 'gap must be >= 0');
  need(Number.isFinite(mass) && mass > 0, 'mass must be > 0');
  need(Number.isFinite(stiffness) && stiffness >= 0, 'stiffness must be >= 0');
  const e = force * gap - 0.5 * stiffness * gap * gap;
  // W(y) = f0 y - k y^2/2 is concave and zero at y = 0, so W(g) > 0 means the section never stops before contact
  const reaches = e > 0;
  return { velocity: reaches ? Math.sqrt(2 * e / mass) : 0, energy: reaches ? e : 0, reaches };
}
