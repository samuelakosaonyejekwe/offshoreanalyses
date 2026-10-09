/**
 * Slug-induced loads on bends, risers, jumpers and topsides piping: momentum force at a bend, force time
 * history of a slug train through a multi-bend spool, dynamic load factor of a support (SDOF), fatigue-cycle
 * list for the fatigue engine, allowable slug velocity and ρv² momentum checks.
 */
import { rk4 } from '../math/ode.js';
import { checkArg } from './driftflux.js';

const area = (D) => Math.PI * D * D / 4;

/**
 * Momentum-flux force of a fluid of density ρ turning through a bend at velocity v:
 * each leg carries ρ A v² along its axis, so the resultant on the bend is
 *   F = ρ A v² √(2 (1 − cosθ))       (θ = 90°: √2 ρ A v²; θ = 180°: 2 ρ A v²),
 * directed along the outward bisector of the bend. The pressure-area term is not included — it is balanced
 * within a continuous pipe and only loads supports at unbalanced ends.
 * Reference: control-volume momentum theorem (e.g. White, Fluid Mechanics, §3.4).
 * @param {object} o
 * @param {number} o.rho fluid density [kg/m³]
 * @param {number} o.v fluid velocity [m/s]
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} [o.bendAngle=π/2] angle the flow is turned through [rad], 0 … π
 * @returns {{axialForce: number, resultant: number, Fx: number, Fy: number}} axialForce = ρ A v² [N];
 *   resultant [N]; components on the bend with x along the incoming flow and y towards the outgoing leg [N]
 *   (Fx = ρ A v² (1 − cosθ), Fy = −ρ A v² sinθ)
 */
export function bendForce({ rho, v, D, bendAngle = Math.PI / 2 }) {
  checkArg('rho', rho, { min: 0 });
  checkArg('v', v);
  checkArg('D', D, { min: 0, exclusiveMin: true });
  checkArg('bendAngle', bendAngle, { min: 0, max: Math.PI });
  const F = rho * area(D) * v * v;
  return { axialForce: F, resultant: F * Math.sqrt(2 * (1 - Math.cos(bendAngle))), Fx: F * (1 - Math.cos(bendAngle)), Fy: -F * Math.sin(bendAngle) };
}

/**
 * Slug load at a bend: the force changes from the film/bubble value to the slug-body value as the slug front
 * arrives, and back as its tail leaves. Slug and film are both taken at the slug (mixture) velocity, the usual
 * conservative assumption, so the step is ΔF = (ρs − ρf) A v² √(2 (1 − cosθ)).
 * Reference: momentum theorem; practice described in e.g. DNV-RP-D101 (2017) §3.12 slugging loads.
 * @param {object} o
 * @param {number} o.rhoSlug slug-body density ρs = HLS ρL + (1 − HLS) ρG [kg/m³]
 * @param {number} o.rhoFilm density of the film/bubble region ρf [kg/m³]
 * @param {number} o.v slug velocity [m/s]
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} [o.bendAngle=π/2] bend angle [rad]
 * @returns {{slugForce: number, filmForce: number, forceRange: number}} resultant with the slug in the bend,
 *   with the film in the bend, and their difference [N]
 */
export function slugBendLoad({ rhoSlug, rhoFilm, v, D, bendAngle = Math.PI / 2 }) {
  checkArg('rhoFilm', rhoFilm, { min: 0 });
  checkArg('rhoSlug', rhoSlug, { min: rhoFilm });
  const s = bendForce({ rho: rhoSlug, v, D, bendAngle }).resultant;
  const f = bendForce({ rho: rhoFilm, v, D, bendAngle }).resultant;
  return { slugForce: s, filmForce: f, forceRange: s - f };
}

/**
 * Force time histories at the bends of a spool for a train of slugs travelling at constant velocity.
 * Each bend sees a trapezoidal pulse per slug: the force follows the fraction of the bend arc occupied by the slug
 * body, so it ramps from the film to the slug value while the front crosses the arc (rise time = arc length / v),
 * holds while the body passes and ramps down at the tail.
 * @param {object} o
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} o.rhoSlug slug-body density [kg/m³]
 * @param {number} o.rhoFilm film-region density [kg/m³]
 * @param {number} o.v slug velocity [m/s], > 0
 * @param {{s: number, bendAngle?: number, radius?: number}[]} o.bends bend positions along the spool s [m],
 *   angles [rad] (default π/2) and centre-line radii [m] (default 1.5 D)
 * @param {{tFront: number, length: number}[]} o.slugs time each slug front passes s = 0 [s] and slug length [m]
 * @param {number} o.duration length of the history [s]
 * @param {number} o.dt sample interval [s]
 * @returns {{t: Float64Array, forces: Float64Array[], peak: number[], range: number[], riseTime: number[]}}
 *   t [s]; resultant force at each bend [N]; peak force, force range (slug − film) [N] and rise time [s] per bend
 */
export function spoolForceHistory({ D, rhoSlug, rhoFilm, v, bends, slugs, duration, dt }) {
  checkArg('v', v, { min: 0, exclusiveMin: true });
  checkArg('duration', duration, { min: 0, exclusiveMin: true });
  checkArg('dt', dt, { min: 0, exclusiveMin: true });
  if (!Array.isArray(bends) || bends.length === 0) throw new RangeError('bends must list at least one bend');
  if (!Array.isArray(slugs)) throw new RangeError('slugs must be an array');
  const n = Math.floor(duration / dt + 1e-9) + 1;
  if (n > 2e6) throw new RangeError('duration/dt gives more than 2e6 samples');
  const t = new Float64Array(n);
  for (let i = 0; i < n; i++) t[i] = i * dt;
  const forces = [];
  const peak = [];
  const range = [];
  const riseTime = [];
  for (const b of bends) {
    checkArg('bend.s', b.s, { min: 0 });
    const ang = b.bendAngle ?? Math.PI / 2;
    const rad = checkArg('bend.radius', b.radius ?? 1.5 * D, { min: 0 });
    const load = slugBendLoad({ rhoSlug, rhoFilm, v, D, bendAngle: ang });
    const arc = rad * ang;
    const tr = arc / v;
    const F = new Float64Array(n).fill(load.filmForce);
    let pk = load.filmForce;
    for (const sl of slugs) {
      checkArg('slug.length', sl.length, { min: 0, exclusiveMin: true });
      const t0 = sl.tFront + b.s / v;
      // Fraction of the bend arc occupied by the slug body: overlap of [tail, front] with the arc [0, La].
      // A slug shorter than the arc never fills it, so its peak load is reduced by Ls/La.
      const i0 = Math.max(0, Math.floor(t0 / dt));
      const i1 = Math.min(n - 1, Math.ceil((t0 + (sl.length + arc) / v) / dt));
      for (let i = i0; i <= i1; i++) {
        const front = v * (t[i] - t0);
        const tail = front - sl.length;
        let w;
        if (arc > 0) w = Math.max(0, Math.min(front, arc) - Math.max(tail, 0)) / arc;
        else w = front >= 0 && tail <= 0 ? 1 : 0;
        const val = load.filmForce + w * load.forceRange;
        if (val > F[i]) F[i] = val;
        if (val > pk) pk = val;
      }
    }
    forces.push(F); peak.push(pk); range.push(load.forceRange); riseTime.push(tr);
  }
  return { t, forces, peak, range, riseTime };
}

/**
 * Dynamic load factor of a support idealised as a single-degree-of-freedom oscillator under a trapezoidal
 * (or rectangular) force pulse: ẍ + 2 ζ ωn ẋ + ωn² x = ωn² f(t), f = pulse normalised to 1, DLF = max |x|.
 * Integrated with classical RK4 over the pulse and `tailPeriods` natural periods of free vibration.
 * Closed-form check (undamped, rectangular pulse of duration td): DLF = 2 sin(π td/Tn) for td/Tn < ½, else 2
 * (Biggs 1964, Introduction to Structural Dynamics, §2.4).
 * @param {object} o
 * @param {number} o.naturalFrequency support natural frequency fn [Hz]
 * @param {number} o.duration plateau duration of the pulse [s] (slug length / velocity)
 * @param {number} [o.riseTime=0] linear rise and decay time [s]
 * @param {number} [o.damping=0] damping ratio ζ [-], 0 … 1
 * @param {number} [o.stepsPerPeriod=200] RK4 steps per natural period (and per rise time) [-]
 * @param {number} [o.tailPeriods=3] free-vibration periods integrated after the pulse [-]
 * @returns {{dlf: number, tPeak: number, staticRatio: number}} dlf = max|x|/x_static [-]; time of the peak
 *   response [s]; staticRatio = td/Tn [-]
 */
export function dynamicLoadFactor({ naturalFrequency, duration, riseTime = 0, damping = 0, stepsPerPeriod = 200, tailPeriods = 3 }) {
  checkArg('naturalFrequency', naturalFrequency, { min: 0, exclusiveMin: true });
  checkArg('duration', duration, { min: 0, exclusiveMin: true });
  checkArg('riseTime', riseTime, { min: 0 });
  checkArg('damping', damping, { min: 0, max: 1 });
  checkArg('stepsPerPeriod', stepsPerPeriod, { min: 20, max: 5000 });
  checkArg('tailPeriods', tailPeriods, { min: 0.5, max: 100 });
  const wn = 2 * Math.PI * naturalFrequency;
  const Tn = 1 / naturalFrequency;
  const tr = riseTime;
  const tPulse = duration + 2 * tr;
  const pulse = (t) => {
    if (tr > 0 && t <= tr) return t / tr;
    if (tr > 0 && t >= tr + duration) return Math.max(0, (tPulse - t) / tr);
    return 1;
  };
  let dlf = 0;
  let tPeak = 0;
  let y = [0, 0];
  let t0 = 0;
  // Integrate each smooth piece of the pulse separately so the kinks fall on step boundaries.
  const pieces = tr > 0 ? [tr, tr + duration, tPulse, tPulse + tailPeriods * Tn] : [duration, duration + tailPeriods * Tn];
  const h = Math.min(Tn, tr > 0 ? tr : Tn) / stepsPerPeriod;
  for (let k = 0; k < pieces.length; k++) {
    const t1 = pieces[k];
    const loaded = k < pieces.length - 1;
    const nSteps = Math.max(1, Math.ceil((t1 - t0) / h));
    if (nSteps > 5e5) throw new RangeError('duration is too long relative to the natural period for time integration (the response is quasi-static: DLF → 1 for slow ramps, 2 for a sudden step)');
    const rhs = (t, s) => [s[1], wn * wn * ((loaded ? pulse(t) : 0) - s[0]) - 2 * damping * wn * s[1]];
    const sol = rk4(rhs, t0, y, t1, nSteps, (t, s) => { const a = Math.abs(s[0]); if (a > dlf) { dlf = a; tPeak = t; } });
    const last = sol.y[sol.y.length - 1];
    y = [last[0], last[1]];
    t0 = t1;
  }
  return { dlf, tPeak, staticRatio: duration / Tn };
}

/**
 * Fatigue-cycle list for slug loading at one bend: each slug gives one force cycle of range DLF × ΔF, converted
 * to a stress range through a user-supplied transfer factor (stress per unit force at the hot spot, from the
 * piping stress model). The slug count per year is frequency × seconds per year × fraction of time in slug
 * flow. With a slug-length distribution the cycles are binned by slug duration (DLF depends on duration).
 * Output format for the fatigue engine: [{ stressRange [Pa], cyclesPerYear }].
 * @param {object} o
 * @param {number} o.forceRange slug − film force range at the bend ΔF [N]
 * @param {number} o.stressPerForce stress at the fatigue hot spot per unit bend force [Pa/N]
 * @param {number} o.slugFrequency slug frequency [1/s]
 * @param {number} [o.slugFlowFraction=1] fraction of the year spent in slug flow [-]
 * @param {number} [o.dlf=1] dynamic load factor applied to every cycle [-] (ignored when `bins` is given)
 * @param {{probability: number, dlf: number}[]} [o.bins] share of slugs and DLF per slug-duration bin [-]
 * @returns {{cycles: {stressRange: number, cyclesPerYear: number}[], cyclesPerYear: number,
 *   maxStressRange: number}} stress ranges [Pa] and annual counts [1/yr]
 */
export function slugFatigueCycles({ forceRange, stressPerForce, slugFrequency, slugFlowFraction = 1, dlf = 1, bins }) {
  checkArg('forceRange', forceRange, { min: 0 });
  checkArg('stressPerForce', stressPerForce, { min: 0 });
  checkArg('slugFrequency', slugFrequency, { min: 0 });
  checkArg('slugFlowFraction', slugFlowFraction, { min: 0, max: 1 });
  const perYear = slugFrequency * 365.25 * 86400 * slugFlowFraction;
  let list;
  if (bins) {
    let sum = 0;
    for (const b of bins) { checkArg('bin.probability', b.probability, { min: 0, max: 1 }); checkArg('bin.dlf', b.dlf, { min: 0 }); sum += b.probability; }
    if (Math.abs(sum - 1) > 1e-6) throw new RangeError(`bins probabilities sum to ${sum}, expected 1`);
    list = bins.map((b) => ({ stressRange: b.dlf * forceRange * stressPerForce, cyclesPerYear: b.probability * perYear }));
  } else {
    checkArg('dlf', dlf, { min: 0 });
    list = [{ stressRange: dlf * forceRange * stressPerForce, cyclesPerYear: perYear }];
  }
  let mx = 0;
  for (const c of list) if (c.stressRange > mx) mx = c.stressRange;
  return { cycles: list, cyclesPerYear: perYear, maxStressRange: mx };
}

/**
 * Allowable slug velocity for an allowable dynamic bend force: inverts DLF (ρs − ρf) A v² √(2 (1 − cosθ)) = Fallow.
 * @param {object} o
 * @param {number} o.allowableForce allowable resultant force range on the bend / support [N]
 * @param {number} o.rhoSlug slug-body density [kg/m³]
 * @param {number} [o.rhoFilm=0] film-region density [kg/m³]
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} [o.bendAngle=π/2] bend angle [rad], > 0
 * @param {number} [o.dlf=2] dynamic load factor [-] (2 = undamped upper bound for a suddenly applied load)
 * @returns {number} allowable slug velocity [m/s]
 */
export function allowableSlugVelocity({ allowableForce, rhoSlug, rhoFilm = 0, D, bendAngle = Math.PI / 2, dlf = 2 }) {
  checkArg('allowableForce', allowableForce, { min: 0 });
  checkArg('rhoFilm', rhoFilm, { min: 0 });
  checkArg('rhoSlug', rhoSlug, { min: rhoFilm, exclusiveMin: true });
  checkArg('D', D, { min: 0, exclusiveMin: true });
  checkArg('bendAngle', bendAngle, { min: 0, max: Math.PI, exclusiveMin: true });
  checkArg('dlf', dlf, { min: 0, exclusiveMin: true });
  return Math.sqrt(allowableForce / (dlf * (rhoSlug - rhoFilm) * area(D) * Math.sqrt(2 * (1 - Math.cos(bendAngle)))));
}

/**
 * Momentum (ρ v²) check for riser bases, separator inlets and inlet devices: ρm vm² against a limit supplied by
 * the user (company or vendor criterion for the nozzle / inlet device in question — no default is assumed).
 * @param {object} o
 * @param {number} o.rho mixture density at the location [kg/m³] (no-slip for inlet nozzles, slug-body density
 *   for slug impact)
 * @param {number} o.v mixture (or slug) velocity [m/s]
 * @param {number} [o.limit] allowable ρ v² [Pa]
 * @returns {{momentum: number, utilisation: number, ok: boolean|null}} momentum = ρ v² [Pa]; utilisation =
 *   momentum/limit [-] and ok (null/NaN when no limit is given)
 */
export function momentumCheck({ rho, v, limit }) {
  checkArg('rho', rho, { min: 0 });
  checkArg('v', v);
  const m = rho * v * v;
  if (limit === undefined) return { momentum: m, utilisation: NaN, ok: null };
  checkArg('limit', limit, { min: 0, exclusiveMin: true });
  return { momentum: m, utilisation: m / limit, ok: m <= limit };
}
