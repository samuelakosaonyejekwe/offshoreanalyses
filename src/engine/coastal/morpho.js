// Bed evolution (Exner equation) and equilibrium-profile relations: Dean profile, Bruun rule, depth of
// closure. SI units.

const G = 9.80665;

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function pos(v, name) {
  need(Number.isFinite(v) && v > 0, `${name} must be a positive finite number`);
}

const LIMITERS = {
  none: () => 1, // unlimited Lax–Wendroff
  minmod: (r) => Math.max(0, Math.min(1, r)),
  vanleer: (r) => (r > 0 ? (2 * r) / (1 + r) : 0),
  superbee: (r) => Math.max(0, Math.min(2 * r, 1), Math.min(r, 2)),
};

/**
 * One-dimensional Exner bed evolution coupled to a user-supplied transport function:
 *   (1 − p) ∂zb/∂t + ∂qs/∂x = 0
 * solved in conservation form on a uniform grid. The face flux is first-order upwind with respect to the
 * bed-form celerity at the face, plus (scheme 'laxWendroff') the flux-limited second-order correction
 * ½ |a| (1 − |a| Δt/Δx) φ(r) Δzb (Sweby 1984). Where the celerity vanishes the upwind direction falls back to
 * the direction of transport.
 *
 * Celerity. The bed-form celerity is c = (∂qs/∂zb)/(1 − p) at a fixed position. With `celerity: 'local'`
 * (default) it is estimated in every cell by one extra call of the transport function with the whole bed raised
 * by a small δ, c_i = [qs_i(zb + δ) − qs_i(zb)] / ((1 − p) δ). The face celerity is the jump ratio
 * Δqs/((1 − p) Δzb) between the two neighbouring cells (the exact shock speed when qs depends on zb only),
 * limited to the range spanned by the two cell values: where the transport also varies with position (changing
 * flow width, depth or forcing) the raw jump ratio is not a celerity at all — it is unbounded where the bed is
 * locally level — and the cell estimate governs. `celerity: 'jump'` uses the raw jump ratio (valid only for
 * qs = f(zb)), and a function (zb, t) → c per cell [m/s] supplies the celerity directly, e.g. when the transport
 * model is invariant to a uniform shift of the bed.
 *
 * Time step (adaptive mode): Δt = cfl·Δx / max|a| over the faces; with a bed-slope term also
 * Δt <= cfl·Δx²(1 − p) / (2 β max|qs|); and never longer than the time in which the fastest-changing cell moves
 * by `maxBedChange`. The last limit is the floor that keeps the step finite when the celerity is zero or tiny
 * while the bed is still changing (transport independent of the bed, supply at a boundary).
 *
 * Bed-slope term (optional): the face flux is reduced going up-slope and increased down-slope,
 * qs → qs − β |qs| ∂zb/∂x, the usual gravitational correction that diffuses bed features (β = `bedSlopeBeta`,
 * typically 1–5; Bailard 1981 gives β = 1/tan φ ≈ 1.6 for bed load).
 * Reference: Exner (1925) Sitzungsber. Akad. Wiss. Wien 134; LeVeque (2002) "Finite Volume Methods for
 * Hyperbolic Problems" §6.11–6.12 and §16 (spatially varying flux); Sweby (1984) SIAM J. Numer. Anal. 21;
 * Bailard (1981) J. Geophys. Res. 86(C11) (slope term).
 * Validity: quasi-steady (decoupled) morphodynamics — the flow is assumed to adapt instantly to the bed, so
 * the transport function receives the current bed and returns qs at cell centres; explicit, stable for
 * Courant number |a| Δt/Δx <= 1. The transport function must not keep a reference to the bed array it is given.
 * @param {object} o
 * @param {number} o.dx cell size [m]
 * @param {ArrayLike<number>} o.zb0 initial bed elevation at cell centres [m]
 * @param {(zb:Float64Array, t:number)=>ArrayLike<number>} o.transport returns the volumetric transport
 *   rate (solid volume) per unit width at each cell centre, positive toward +x [m²/s]
 * @param {number} o.tEnd end time [s]
 * @param {number} [o.porosity] bed porosity p [-], default 0.4
 * @param {'upwind'|'laxWendroff'} [o.scheme] default 'laxWendroff'
 * @param {'none'|'minmod'|'vanleer'|'superbee'} [o.limiter] default 'vanleer'
 * @param {number} [o.cfl] Courant number for the adaptive step [-], default 0.8
 * @param {number} [o.dt] fixed time step [s] (overrides cfl and maxBedChange; the last step is shortened to hit
 *   tEnd)
 * @param {number} [o.dtMax] upper limit on the adaptive step [s]
 * @param {'local'|'jump'|((zb:Float64Array, t:number)=>ArrayLike<number>)} [o.celerity] celerity estimate,
 *   default 'local' (see above); a function returns the celerity per cell [m/s]
 * @param {number} [o.maxBedChange] largest bed-level change of any cell in one adaptive step [m]; default
 *   cfl × max(initial bed relief, 0.001 Δx)
 * @param {number} [o.bedSlopeBeta] bed-slope coefficient β [-], default 0 (no slope term)
 * @param {'periodic'|'zeroGradient'|{q:number|((t:number)=>number)}} [o.left] upstream boundary: periodic,
 *   zero-gradient (flux equals the edge-cell transport) or prescribed sediment flux [m²/s]
 * @param {'periodic'|'zeroGradient'|{q:number|((t:number)=>number)}} [o.right] as left
 * @param {number} [o.maxSteps] default 1e7
 * @returns {{zb:Float64Array,t:number,steps:number,courantMax:number,dtMin:number,dtMax:number,
 *   mass:{initial:number,final:number,boundaryNet:number,error:number,relativeError:number}}}
 *   final bed [m]; courantMax = largest |a| Δt/Δx used [-]; dtMin, dtMax = smallest and largest step taken [s];
 *   mass = solid volume per unit width above datum (1 − p) Σ zb Δx [m²], boundaryNet =
 *   time-integrated inflow minus outflow [m²], error = final − initial − boundaryNet
 */
export function exnerEvolve({
  dx, zb0, transport, tEnd, porosity = 0.4, scheme = 'laxWendroff', limiter = 'vanleer', cfl = 0.8, dt, dtMax,
  left = 'zeroGradient', right = 'zeroGradient', maxSteps = 1e7, celerity = 'local', maxBedChange, bedSlopeBeta = 0,
}) {
  pos(dx, 'dx');
  need(zb0 && zb0.length >= 3, 'zb0 must have at least 3 cells');
  need(typeof transport === 'function', 'transport must be a function');
  need(Number.isFinite(tEnd) && tEnd >= 0, 'tEnd must be >= 0');
  need(porosity >= 0 && porosity < 1, 'porosity must lie in [0, 1)');
  need(scheme === 'upwind' || scheme === 'laxWendroff', "scheme must be 'upwind' or 'laxWendroff'");
  const phi = LIMITERS[limiter];
  need(phi !== undefined, 'limiter must be none, minmod, vanleer or superbee');
  need(cfl > 0 && cfl <= 1, 'cfl must lie in (0, 1]');
  if (dt !== undefined) pos(dt, 'dt');
  if (dtMax !== undefined) pos(dtMax, 'dtMax');
  need(celerity === 'local' || celerity === 'jump' || typeof celerity === 'function', "celerity must be 'local', 'jump' or a function");
  if (maxBedChange !== undefined) pos(maxBedChange, 'maxBedChange');
  need(Number.isFinite(bedSlopeBeta) && bedSlopeBeta >= 0, 'bedSlopeBeta must be >= 0');
  const periodic = left === 'periodic' || right === 'periodic';
  need(!periodic || (left === 'periodic' && right === 'periodic'), 'periodic boundaries must be set on both sides');
  const bFlux = (b, name) => {
    if (b === 'periodic' || b === 'zeroGradient') return null;
    need(b && b.q !== undefined, `${name} must be 'periodic', 'zeroGradient' or { q }`);
    return typeof b.q === 'function' ? b.q : () => b.q;
  };
  const qLeft = bFlux(left, 'left');
  const qRight = bFlux(right, 'right');
  const n = zb0.length;
  const z = Float64Array.from(zb0);
  const F = new Float64Array(n + 1); // face fluxes of solid volume [m²/s]; face i lies between cells i−1 and i
  const a = new Float64Array(n + 1);
  const dz = new Float64Array(n + 1);
  const cc = new Float64Array(n); // cell celerity [m/s]
  const zPert = new Float64Array(n);
  const pf = 1 / (1 - porosity);
  const second = scheme === 'laxWendroff';
  const useCell = celerity !== 'jump';
  let zLo = Infinity;
  let zHi = -Infinity;
  let zAbs = 0;
  for (let i = 0; i < n; i++) {
    need(Number.isFinite(z[i]), 'zb0 must be finite');
    if (z[i] < zLo) zLo = z[i];
    if (z[i] > zHi) zHi = z[i];
    if (Math.abs(z[i]) > zAbs) zAbs = Math.abs(z[i]);
  }
  const dzLimit = maxBedChange !== undefined ? maxBedChange : cfl * Math.max(zHi - zLo, 1e-3 * dx);
  // Perturbation for the celerity estimate: small against any depth scale, large against round-off.
  const delta = 1e-6 * Math.max(zAbs, 1);
  const mass = () => {
    let s = 0;
    for (let i = 0; i < n; i++) s += z[i];
    return s * dx * (1 - porosity);
  };
  const m0 = mass();
  let boundaryNet = 0;
  let t = 0;
  let steps = 0;
  let courantMax = 0;
  let dtLo = Infinity;
  let dtHi = 0;
  const lo = periodic ? 0 : 1;
  const hi = n - 1;
  const fluxes = (q, step) => {
    for (let f = lo; f <= hi; f++) {
      const il = f === 0 ? n - 1 : f - 1;
      const af = a[f];
      const dir = af !== 0 ? af : q[f] + q[il];
      let flux = dir >= 0 ? q[il] : q[f];
      if (second && af !== 0) {
        // Ratio of the upwind bed jump to the local one.
        let fu = dir >= 0 ? f - 1 : f + 1;
        let r = 0;
        let ok = dz[f] !== 0;
        if (periodic) fu = (fu + n) % n;
        else if (fu < 1 || fu > n - 1) ok = false;
        if (ok) r = dz[fu] / dz[f];
        const lim = ok ? phi(r) : 0;
        const aa = Math.abs(af);
        flux += 0.5 * aa * (1 - (aa * step) / dx) * lim * dz[f] / pf;
      }
      if (bedSlopeBeta > 0) flux -= bedSlopeBeta * 0.5 * (Math.abs(q[f]) + Math.abs(q[il])) * dz[f] / dx;
      F[f] = flux;
    }
    if (periodic) F[n] = F[0];
    else {
      F[0] = qLeft ? qLeft(t) : q[0];
      F[n] = qRight ? qRight(t) : q[n - 1];
    }
  };
  while (t < tEnd) {
    need(steps < maxSteps, 'maxSteps exceeded');
    if (celerity === 'local') {
      // the transport function may reuse its output buffer, so the perturbed result is consumed before the
      // unperturbed call
      for (let i = 0; i < n; i++) zPert[i] = z[i] + delta;
      const qp = transport(zPert, t);
      need(qp && qp.length === n, 'transport must return one value per cell');
      for (let i = 0; i < n; i++) cc[i] = qp[i];
    } else if (useCell) {
      const cu = celerity(z, t);
      need(cu && cu.length === n, 'celerity must return one value per cell');
      for (let i = 0; i < n; i++) cc[i] = cu[i];
    }
    const q = transport(z, t);
    need(q && q.length === n, 'transport must return one value per cell');
    let amax = 0;
    let qmax = 0;
    for (let i = 0; i < n; i++) {
      need(Number.isFinite(q[i]), 'transport returned a non-finite value');
      if (Math.abs(q[i]) > qmax) qmax = Math.abs(q[i]);
      if (celerity === 'local') {
        // the difference quotient carries round-off of order 1e-9; rounding to 7 significant figures removes it,
        // so that a uniform celerity gives a uniform Courant number
        const c = (pf * (cc[i] - q[i])) / delta;
        cc[i] = c === 0 || !Number.isFinite(c) ? c : Number(c.toPrecision(7));
      }
      if (useCell) need(Number.isFinite(cc[i]), 'celerity is not finite');
    }
    for (let f = lo; f <= hi; f++) {
      const il = f === 0 ? n - 1 : f - 1;
      const d = z[f] - z[il];
      const dq = q[f] - q[il];
      dz[f] = d;
      const scale = 1e-12 * (Math.abs(z[f]) + Math.abs(z[il]) + 1e-30);
      const level = !(Math.abs(d) > scale);
      let af = level ? 0 : (pf * dq) / d;
      if (useCell) {
        const cA = cc[il];
        const cB = cc[f];
        if (level) af = 0.5 * (cA + cB);
        else {
          const w = 1e-6 * (Math.abs(cA) + Math.abs(cB));
          const cLo = Math.min(cA, cB) - w;
          const cHi = Math.max(cA, cB) + w;
          af = af < cLo ? cLo : af > cHi ? cHi : af;
        }
      }
      a[f] = af;
      if (Math.abs(af) > amax) amax = Math.abs(af);
    }
    let step = dt !== undefined ? dt : amax > 0 ? (cfl * dx) / amax : tEnd - t;
    if (dt === undefined && bedSlopeBeta > 0 && qmax > 0) {
      const sd = (cfl * dx * dx) / (2 * pf * bedSlopeBeta * qmax);
      if (step > sd) step = sd;
    }
    if (dtMax !== undefined && step > dtMax) step = dtMax;
    if (t + step > tEnd) step = tEnd - t;
    fluxes(q, step);
    if (dt === undefined) {
      let rate = 0;
      for (let i = 0; i < n; i++) {
        const r = Math.abs(F[i + 1] - F[i]);
        if (r > rate) rate = r;
      }
      rate *= pf / dx; // fastest bed-level change [m/s]
      if (rate * step > dzLimit) {
        step = dzLimit / rate;
        fluxes(q, step);
      }
    }
    let hit = false;
    if (t + step >= tEnd - 1e-12 * Math.max(tEnd, 1)) {
      step = tEnd - t;
      hit = true;
    }
    const nu = (amax * step) / dx;
    if (nu > courantMax) courantMax = nu;
    const c = (step * pf) / dx;
    for (let i = 0; i < n; i++) z[i] -= c * (F[i + 1] - F[i]);
    boundaryNet += step * (F[0] - F[n]);
    t = hit ? tEnd : t + step;
    steps++;
    if (step < dtLo) dtLo = step;
    if (step > dtHi) dtHi = step;
  }
  const mf = mass();
  const error = mf - m0 - boundaryNet;
  return {
    zb: z, t, steps, courantMax, dtMin: steps > 0 ? dtLo : 0, dtMax: dtHi,
    mass: { initial: m0, final: mf, boundaryNet, error, relativeError: Math.abs(m0) > 0 ? error / Math.abs(m0) : error },
  };
}

/**
 * Dean equilibrium profile scale parameter from the sediment fall velocity:
 * A = 0.067 ws^0.44 with ws in cm/s and A in m^(1/3).
 * Reference: Dean (1987) "Coastal sediment processes: toward engineering solutions", Proc. Coastal
 * Sediments '87 (fit to Moore's 1982 data). Validity: sand, ws about 1–10 cm/s.
 * @param {number} ws fall velocity [m/s]
 * @returns {number} A [m^(1/3)]
 */
export function deanParameter(ws) {
  pos(ws, 'ws');
  return 0.067 * Math.pow(100 * ws, 0.44);
}

/**
 * Dean equilibrium beach profile h(y) = A y^(2/3) (uniform wave-energy dissipation per unit volume).
 * Reference: Bruun (1954); Dean (1977) Ocean Eng. Rep. 12, Univ. Delaware; Dean (1991) J. Coastal Res. 7(1).
 * Validity: sandy profile seaward of the shoreline out to the depth of closure; infinite slope at y = 0.
 * @param {object} o
 * @param {number|ArrayLike<number>} o.y distance(s) offshore from the shoreline [m]
 * @param {number} o.A profile scale parameter [m^(1/3)] (see deanParameter)
 * @returns {{h:number|Float64Array, widthToDepth:(h:number)=>number}} depth(s) [m]; widthToDepth(h) =
 *   (h/A)^(3/2) [m]
 */
export function deanProfile({ y, A }) {
  pos(A, 'A');
  const one = (v) => {
    need(Number.isFinite(v) && v >= 0, 'y must be >= 0');
    return A * Math.pow(v, 2 / 3);
  };
  const h = typeof y === 'number' ? one(y) : Float64Array.from(y, one);
  return { h, widthToDepth: (d) => Math.pow(d / A, 1.5) };
}

/**
 * Bruun rule: shoreline retreat R = S L /(B + hc) for a sea-level rise S, with L the cross-shore width of
 * the active profile, B the berm (or dune) height and hc (= h*) the depth of closure. If L is omitted it is
 * taken from a Dean profile, L = (hc/A)^(3/2) with hc = h*.
 * Reference: Bruun (1962) J. Waterways Harbors Div. ASCE 88(WW1). Validity: two-dimensional equilibrium
 * profile translation with a closed sediment budget — no longshore gradients, no overwash, no rock.
 * @param {object} o
 * @param {number} o.S sea-level rise [m]
 * @param {number} o.hClosure depth of closure h* [m]
 * @param {number} o.B berm height above sea level [m]
 * @param {number} [o.L] active profile width [m]
 * @param {number} [o.A] Dean parameter used when L is omitted [m^(1/3)]
 * @returns {{retreat:number, L:number, slope:number, ratio:number}} retreat [m], L [m], mean active-profile
 *   slope (B + h*)/L [-], retreat per unit rise [-]
 */
export function bruunRule({ S, hClosure, B, L, A }) {
  need(Number.isFinite(S) && S >= 0, 'S must be >= 0');
  pos(hClosure, 'hClosure');
  need(Number.isFinite(B) && B >= 0, 'B must be >= 0');
  let width = L;
  if (width === undefined) {
    pos(A, 'A (required when L is omitted)');
    width = Math.pow(hClosure / A, 1.5);
  }
  pos(width, 'L');
  const slope = (B + hClosure) / width;
  return { retreat: S / slope, L: width, slope, ratio: 1 / slope };
}

/**
 * Depth of closure.
 *   'hallermeier' (inner limit): hc = 2.28 He − 68.5 He²/(g Te²)
 *   'birkemeier':                hc = 1.75 He − 57.9 He²/(g Te²)
 * with He the nearshore significant wave height exceeded 12 hours per year and Te its period.
 * Reference: Hallermeier (1978) Proc. 16th ICCE; Hallermeier (1981) Coastal Eng. 4; Birkemeier (1985)
 * J. Waterway Port Coastal Ocean Eng. 111(3). Validity: quartz sand, open wave-dominated coasts.
 * @param {object} o
 * @param {number} o.He extreme significant wave height (12 h/yr) [m]
 * @param {number} o.Te associated period [s]
 * @param {'hallermeier'|'birkemeier'} [o.method] default 'hallermeier'
 * @param {number} [o.g] gravity [m/s²]
 * @returns {number} depth of closure below mean low water [m]
 */
export function depthOfClosure({ He, Te, method = 'hallermeier', g = G }) {
  pos(He, 'He');
  pos(Te, 'Te');
  const s = (He * He) / (g * Te * Te);
  if (method === 'hallermeier') return 2.28 * He - 68.5 * s;
  if (method === 'birkemeier') return 1.75 * He - 57.9 * s;
  throw new RangeError("method must be 'hallermeier' or 'birkemeier'");
}

/**
 * Hallermeier (1981) outer (seaward) limit of the shoal zone: hi = (Hs,mean − 0.3 σs) Ts,mean (g/(5000 d50))^0.5.
 * Reference: Hallermeier (1981) Coastal Eng. 4, 253–277. Validity: quartz sand 0.16–0.42 mm.
 * @param {object} o
 * @param {number} o.HsMean annual mean significant wave height [m]
 * @param {number} o.sigmaHs standard deviation of Hs [m]
 * @param {number} o.TsMean annual mean significant period [s]
 * @param {number} o.d50 median grain size [m]
 * @param {number} [o.g] gravity [m/s²]
 * @returns {number} outer closure depth [m]
 */
export function depthOfClosureOuter({ HsMean, sigmaHs, TsMean, d50, g = G }) {
  pos(HsMean, 'HsMean');
  need(Number.isFinite(sigmaHs) && sigmaHs >= 0, 'sigmaHs must be >= 0');
  pos(TsMean, 'TsMean');
  pos(d50, 'd50');
  return (HsMean - 0.3 * sigmaHs) * TsMean * Math.sqrt(g / (5000 * d50));
}
