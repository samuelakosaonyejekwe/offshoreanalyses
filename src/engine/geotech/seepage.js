/**
 * Seepage: one-dimensional layered Darcy flow, two-dimensional confined seepage under structures with
 * cut-off walls (finite volumes), piping and heave safety, and wave-induced pore pressure in the seabed.
 *
 * SI units: m, s, m/s (hydraulic conductivity), Pa, N/m^3; heads in m.
 *
 * References
 *   Darcy H. (1856); Terzaghi K. (1943) Theoretical Soil Mechanics, ch. 12 (piping, heave prism).
 *   Harr M.E. (1962) Groundwater and Seepage, McGraw-Hill (sheet-pile and weir solutions).
 *   Putnam J.A. (1949) Loss of wave energy due to percolation in a permeable sea bottom, Trans. AGU 30(3).
 *   Yamamoto T., Koning H.L., Sellmeijer H. & van Hijum E. (1978) J. Fluid Mech. 87; Madsen O.S. (1978)
 *     Geotechnique 28(4) – for a deep, stiff, saturated bed both reduce to p = p0 e^(−kz).
 *   MacPherson H. (1978) Wave forces on pipeline buried in permeable seabed, JWPCOE 104(WW4).
 *   Burcharth H.F., Liu Z. & Troch P. (1999) Scaling of core material in rubble mound breakwater model
 *     tests, COPEDEC V.
 */
import { CooBuilder, pcg } from '../math/sparse.js';

function need(name, v, min = -Infinity, max = Infinity) {
  if (typeof v !== 'number' || Number.isNaN(v) || v < min || v > max) {
    throw new RangeError(`${name} must be a number in [${min}, ${max}] (got ${v})`);
  }
  return v;
}

/**
 * Steady one-dimensional Darcy flow through layers in series:
 *   v = Δh / Σ (Li/ki),  equivalent conductivity normal to the layering k⊥ = ΣLi / Σ(Li/ki),
 *   parallel to the layering k∥ = Σ(ki Li)/ΣLi; head loss in layer i = v Li/ki.
 * @param {object} o
 * @param {{thickness:number,k:number}[]} o.layers in flow order; thickness [m], k [m/s]
 * @param {number} o.headDrop total head difference [m]
 * @param {number} [o.area=1] flow area [m^2]
 * @returns {{velocity:number,Q:number,kNormal:number,kParallel:number,headLoss:number[],gradients:number[],
 *   heads:number[]}} Darcy velocity [m/s]; discharge [m^3/s]; head loss and gradient per layer; head
 *   (relative to the outflow face) at each interface starting with the inflow face [m]
 */
export function layeredDarcy({ layers, headDrop, area = 1 }) {
  if (!Array.isArray(layers) || layers.length === 0) throw new RangeError('layers must be a non-empty array');
  need('area', area, 0);
  let sumL = 0;
  let sumLk = 0;
  let sumKL = 0;
  for (const l of layers) {
    need('layer.thickness', l.thickness, 1e-12);
    need('layer.k', l.k, 1e-30);
    sumL += l.thickness;
    sumLk += l.thickness / l.k;
    sumKL += l.thickness * l.k;
  }
  const velocity = headDrop / sumLk;
  const headLoss = layers.map((l) => (velocity * l.thickness) / l.k);
  const heads = [headDrop];
  for (const dh of headLoss) heads.push(heads[heads.length - 1] - dh);
  return { velocity, Q: velocity * area, kNormal: sumL / sumLk, kParallel: sumKL / sumL, headLoss, gradients: headLoss.map((dh, i) => dh / layers[i].thickness), heads };
}

/**
 * Critical hydraulic gradient for upward flow (quick condition): i_cr = γ′/γw = (Gs − 1)/(1 + e).
 * @param {object} o gammaSub [N/m^3], gammaW [N/m^3] (default 9810)
 * @returns {number} i_cr [-]
 */
export function criticalGradient({ gammaSub, gammaW = 9810 }) {
  need('gammaSub', gammaSub, 0);
  need('gammaW', gammaW, 1e-9);
  return gammaSub / gammaW;
}

/**
 * Two-dimensional steady confined seepage in a rectangular soil layer below a structure and/or cut-off
 * walls: ∂/∂x(kx ∂h/∂x) + ∂/∂z(kz ∂h/∂z) = 0, solved on a cell-centred finite-volume grid (second order)
 * with a preconditioned conjugate-gradient solver.
 * Geometry (x to the right, z downward from the soil surface): soil occupies 0 ≤ x ≤ width,
 * 0 ≤ z ≤ thickness over an impervious base; the lateral boundaries are impervious (place them several
 * layer thicknesses away). On the surface the head is `headUp` for x < structure.x1, `headDown` for
 * x > structure.x2, and the surface is impervious in between (base of the structure). Cut-off walls
 * (sheet piles) are impervious vertical lines from the surface to depth `depth` at position `x`.
 * Without a structure give x1 = x2 = position of the sheet pile.
 * Outputs: head field, discharge per metre run, uplift pressure under the structure, maximum exit
 * gradient on the downstream surface, and factors of safety against piping (i_cr / i_exit) and against
 * heave of Terzaghi's prism (depth d × width d/2 on the downstream side of the last cut-off):
 * FS = γ′ d / (γw h_avg), h_avg the mean excess head at the wall-tip level across the prism.
 * The exit gradient is evaluated over the top half-cell and converges slowly near a wall singularity:
 * refine the grid until it stabilises.
 * @param {object} o
 * @param {number} o.width [m]  @param {number} o.thickness [m]
 * @param {number} [o.nx=120]  @param {number} [o.nz=40] cells
 * @param {number} o.kx horizontal conductivity [m/s]  @param {number} [o.kz=o.kx] vertical [m/s]
 * @param {number} o.headUp upstream head above the soil surface [m]
 * @param {number} [o.headDown=0] downstream head [m]
 * @param {{x1:number,x2:number}} o.structure extent of the impervious base [m]
 * @param {{x:number,depth:number}[]} [o.cutoffs=[]] sheet piles [m]
 * @param {number} [o.gammaW=9810] [N/m^3]  @param {number} [o.gammaSub] γ′ for the safety factors [N/m^3]
 * @returns {{head:Float64Array,nx:number,nz:number,dx:number,dz:number,q:number,qOut:number,
 *   exitGradient:number,xExit:number,uplift:{x:number[],pressure:number[],force:number,xBar:number},
 *   fsPiping:number|undefined,fsHeave:number|undefined,converged:boolean}}
 *   head[j*nx + i] at cell centres (row j from the surface) [m]; q inflow and qOut outflow per metre
 *   [m^2/s]; uplift pressure [Pa] under the structure relative to the base level, force [N/m] and its
 *   position [m]
 */
export function confinedSeepage2D(o) {
  const { width, thickness, nx = 120, nz = 40, kx, headUp, headDown = 0, structure, cutoffs = [], gammaW = 9810, gammaSub } = o;
  need('width', width, 1e-9);
  need('thickness', thickness, 1e-9);
  need('kx', kx, 1e-30);
  const kz = o.kz === undefined ? kx : need('kz', o.kz, 1e-30);
  if (!Number.isInteger(nx) || !Number.isInteger(nz) || nx < 4 || nz < 2) throw new RangeError('nx ≥ 4 and nz ≥ 2 must be integers');
  if (!structure) throw new RangeError('structure {x1, x2} is required');
  need('structure.x1', structure.x1, 0, width);
  need('structure.x2', structure.x2, structure.x1, width);
  const dx = width / nx;
  const dz = thickness / nz;
  const n = nx * nz;
  // blocked vertical faces: face f (between columns f−1 and f) is blocked for rows j < rows[f]
  const blockedRows = new Int32Array(nx + 1);
  for (const c of cutoffs) {
    need('cutoff.x', c.x, 0, width);
    need('cutoff.depth', c.depth, 0, thickness * (1 - 1e-9));
    const f = Math.min(Math.max(Math.round(c.x / dx), 1), nx - 1);
    blockedRows[f] = Math.max(blockedRows[f], Math.round(c.depth / dz));
  }
  const cxF = (kx * dz) / dx; // conductance between horizontal neighbours
  const czF = (kz * dx) / dz;
  const cTop = (kz * dx) / (dz / 2);
  const topType = new Int8Array(nx); // 1 upstream, 2 downstream, 0 impervious
  for (let i = 0; i < nx; i++) {
    const xm = (i + 0.5) * dx;
    topType[i] = xm < structure.x1 ? 1 : xm > structure.x2 ? 2 : 0;
  }
  if (!topType.includes(1) || !topType.includes(2)) throw new RangeError('the structure must leave open soil surface on both sides');
  const coo = new CooBuilder(n);
  const b = new Float64Array(n);
  for (let j = 0; j < nz; j++) {
    for (let i = 0; i < nx; i++) {
      const id = j * nx + i;
      let diag = 0;
      if (i > 0 && j >= blockedRows[i]) { coo.add(id, id - 1, -cxF); diag += cxF; }
      if (i < nx - 1 && j >= blockedRows[i + 1]) { coo.add(id, id + 1, -cxF); diag += cxF; }
      if (j > 0) { coo.add(id, id - nx, -czF); diag += czF; }
      if (j < nz - 1) { coo.add(id, id + nx, -czF); diag += czF; }
      if (j === 0 && topType[i] !== 0) { diag += cTop; b[id] = cTop * (topType[i] === 1 ? headUp : headDown); }
      coo.add(id, id, diag);
    }
  }
  const x0 = new Float64Array(n).fill(0.5 * (headUp + headDown));
  const sol = pcg(coo.toCsr(), b, { x0, tol: 1e-11, maxIter: 20 * n, precond: 'ic0' });
  const head = sol.x;
  let q = 0;
  let qOut = 0;
  let exitGradient = 0;
  let xExit = 0;
  const ux = [];
  const up = [];
  let force = 0;
  let mom = 0;
  for (let i = 0; i < nx; i++) {
    const xm = (i + 0.5) * dx;
    if (topType[i] === 1) q += cTop * (headUp - head[i]);
    else if (topType[i] === 2) {
      qOut += cTop * (head[i] - headDown);
      const g = (head[i] - headDown) / (dz / 2);
      if (g > exitGradient) { exitGradient = g; xExit = xm; }
    } else {
      // zero-flux surface: the head at the base of the structure equals the cell value to second order
      const hTop = nz > 1 ? head[i] + (head[i] - head[nx + i]) / 8 : head[i];
      const p = gammaW * hTop;
      ux.push(xm); up.push(p); force += p * dx; mom += p * dx * xm;
    }
  }
  let fsPiping;
  let fsHeave;
  if (gammaSub !== undefined) {
    need('gammaSub', gammaSub, 0);
    fsPiping = exitGradient > 0 ? gammaSub / gammaW / exitGradient : Infinity;
    if (cutoffs.length > 0) {
      let last = cutoffs[0];
      for (const c of cutoffs) if (c.x > last.x) last = c;
      const f = Math.min(Math.max(Math.round(last.x / dx), 1), nx - 1);
      const rows = blockedRows[f];
      if (rows > 0) {
        const d = rows * dz;
        const cols = Math.max(Math.round(d / 2 / dx), 1);
        let s = 0;
        let cnt = 0;
        // head at the tip level = interface between rows (rows − 1) and rows
        for (let i = f; i < Math.min(f + cols, nx); i++) {
          const hTip = rows < nz ? 0.5 * (head[(rows - 1) * nx + i] + head[rows * nx + i]) : head[(rows - 1) * nx + i];
          s += hTip - headDown;
          cnt++;
        }
        const hAvg = s / cnt;
        fsHeave = hAvg > 0 ? (gammaSub * d) / (gammaW * hAvg) : Infinity;
      }
    }
  }
  return {
    head, nx, nz, dx, dz, q, qOut, exitGradient, xExit,
    uplift: { x: ux, pressure: up, force, xBar: force > 0 ? mom / force : NaN },
    fsPiping, fsHeave, converged: sol.converged,
  };
}

/**
 * Wave number from the linear dispersion relation ω² = g k tanh(k d) (Newton iteration).
 * @param {number} T wave period [s]
 * @param {number} d water depth [m]
 * @returns {number} k [rad/m]
 */
export function waveNumber(T, d) {
  need('T', T, 1e-9);
  need('d', d, 1e-9);
  const g = 9.80665;
  const w2 = (2 * Math.PI / T) ** 2;
  let k = w2 / g / Math.sqrt(Math.tanh((w2 * d) / g));
  for (let i = 0; i < 50; i++) {
    const th = Math.tanh(k * d);
    const f = g * k * th - w2;
    const df = g * th + g * k * d * (1 - th * th);
    const dk = f / df;
    k -= dk;
    if (Math.abs(dk) < 1e-14 * k) break;
  }
  return k;
}

/**
 * Wave-induced pore-pressure amplitude in a porous seabed (linear waves).
 *   Bed pressure amplitude:  p0 = γw H / (2 cosh k d).
 *   Infinitely thick, isotropic, fully saturated bed with a stiff skeleton (Yamamoto et al. 1978; Madsen
 *   1978 limit – the pore pressure is then independent of permeability):  p(z) = p0 e^(−k z).
 *   Finite thickness h over an impervious rigid base (Putnam 1949, incompressible potential flow):
 *   p(z) = p0 cosh(k (h − z)) / cosh(k h).
 * Partly saturated or soft beds show additional attenuation and phase lag not modelled here.
 * @param {object} o
 * @param {number} o.H wave height [m]  @param {number} o.d water depth [m]
 * @param {number} [o.T] period [s]  @param {number} [o.k] wave number [rad/m] (one of T, k)
 * @param {number} o.z depth below the seabed [m]
 * @param {number} [o.thickness=Infinity] bed thickness [m]
 * @param {number} [o.gammaW=10055] [N/m^3]
 * @returns {{p0:number,p:number,dpdz:number,k:number,ratio:number}} amplitudes: bed pressure [Pa], pore
 *   pressure at z [Pa], vertical pressure-gradient magnitude at z [Pa/m]; k [rad/m]; p/p0 [-]
 */
export function wavePorePressure({ H, d, T, k, z, thickness = Infinity, gammaW = 10055 }) {
  need('H', H, 0);
  need('d', d, 1e-9);
  need('z', z, 0);
  const kk = k !== undefined ? need('k', k, 1e-12) : waveNumber(need('T', T, 1e-9), d);
  const p0 = (gammaW * H) / (2 * Math.cosh(kk * d));
  let ratio;
  let grad;
  if (Number.isFinite(thickness)) {
    need('thickness', thickness, z);
    ratio = Math.cosh(kk * (thickness - z)) / Math.cosh(kk * thickness);
    grad = (kk * Math.sinh(kk * (thickness - z))) / Math.cosh(kk * thickness);
  } else {
    ratio = Math.exp(-kk * z);
    grad = kk * ratio;
  }
  return { p0, p: p0 * ratio, dpdz: p0 * grad, k: kk, ratio };
}

/**
 * Wave-induced seepage (uplift) force amplitude on a buried pipeline, per unit length. For a pipe small
 * compared with the wavelength the undisturbed pore-pressure field around it is a uniform gradient; an
 * impervious circular inclusion in a potential field doubles the pressure variation on its surface
 * (φ = −G (r + a²/r) cos θ), so
 *     F = inclusionFactor · (π D²/4) · |∂p/∂z|,   inclusionFactor = 2 (MacPherson 1978),
 * acting upward under the wave trough. Compare with the submerged weight plus soil resistance.
 * @param {object} o arguments of wavePorePressure with z = depth of the pipe axis, plus
 * @param {number} o.D pipe outer diameter [m]
 * @param {number} [o.inclusionFactor=2] [-]
 * @returns {{F:number,dpdz:number,p0:number}} uplift amplitude [N/m]; gradient [Pa/m]; bed pressure [Pa]
 */
export function buriedPipeWaveUplift(o) {
  const { D, inclusionFactor = 2 } = o;
  need('D', D, 1e-9);
  need('inclusionFactor', inclusionFactor, 0);
  const r = wavePorePressure(o);
  return { F: (inclusionFactor * Math.PI * D * D * r.dpdz) / 4, dpdz: r.dpdz, p0: r.p0 };
}

/**
 * Attenuation of the wave-induced pore-pressure amplitude inside the core of a rubble-mound breakwater,
 *     p(x) = p0 exp(−δ · 2π x / L′),
 * x being the horizontal distance into the core. The damping coefficient δ is an argument. Burcharth,
 * Liu & Troch (1999) proposed δ = a0 n^(1/2) Lp² / (Hs b) with a0 ≈ 0.014 and L′ = Lp / √D′ (D′ ≈ 1.4 for
 * quarry rock, 1.5 for sand); `burcharthDelta` evaluates that expression with a0 as an argument – check
 * the coefficient against the source before design use.
 * @param {object} o
 * @param {number} o.p0 pressure amplitude at the core boundary [Pa]
 * @param {number} o.x distance into the core [m]
 * @param {number} o.delta damping coefficient [-]
 * @param {number} o.Lprime wavelength in the core [m]
 * @returns {number} p(x) [Pa]
 */
export function coreAttenuation({ p0, x, delta, Lprime }) {
  need('x', x, 0);
  need('delta', delta, 0);
  need('Lprime', Lprime, 1e-9);
  return p0 * Math.exp((-delta * 2 * Math.PI * x) / Lprime);
}

/**
 * Damping coefficient of Burcharth et al. (1999): δ = a0 √n Lp² / (Hs b), and L′ = Lp / √D′.
 * @param {object} o
 * @param {number} o.n core porosity [-]  @param {number} o.Lp peak wavelength [m]
 * @param {number} o.Hs significant wave height [m]  @param {number} o.b core width at the level considered [m]
 * @param {number} [o.a0=0.014] empirical coefficient [-]  @param {number} [o.Dprime=1.4] [-]
 * @returns {{delta:number,Lprime:number}}
 */
export function burcharthDelta({ n, Lp, Hs, b, a0 = 0.014, Dprime = 1.4 }) {
  need('n', n, 0, 1);
  need('Lp', Lp, 1e-9);
  need('Hs', Hs, 1e-9);
  need('b', b, 1e-9);
  return { delta: (a0 * Math.sqrt(n) * Lp * Lp) / (Hs * b), Lprime: Lp / Math.sqrt(Dprime) };
}
