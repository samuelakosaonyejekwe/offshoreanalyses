// One-dimensional axial pipe–soil model for walking: an elastic bar (EA) on elastic–perfectly-plastic friction
// sliders, loaded by a sequence of free-strain fields and cycled. Same idealisation as Carr, Sinclair & Bruton (2006),
// OTC 17945, whose closed forms are its rigid-plastic limits.
//
// Each load step minimises the (convex, piecewise-quadratic) potential energy. The solver here is a Newton iteration
// with an exact line search on the directional derivative, which is monotone for a convex energy, so every iteration
// lowers the energy and the active set of sliding nodes cannot cycle; a hard iteration cap ends a step that will not
// converge with an error instead of subdividing it.
import { solveTridiagonal } from '../../engine/math/linalg.js';

/**
 * @param {object} o
 * @param {number} o.L length [m] @param {number} o.EA axial stiffness [N] @param {number} o.W submerged weight [N/m]
 * @param {number} o.mu axial friction coefficient [-]
 * @param {Float64Array[]} o.steps free axial strain per element for each load step of one cycle
 * @param {number} [o.cycles=3] @param {number} [o.slopeRad=0] seabed slope, positive when it falls towards +x [rad]
 * @param {number} [o.tensionAtEnd=0] constant tension at x = L pulling towards +x [N]
 * @param {number} [o.tensionAtStart=0] constant tension at x = 0 pulling towards −x [N]
 * @param {number} [o.uMob=0.002] friction mobilisation displacement [m]
 * @param {{x:number, stiffness:number}} [o.anchor] linear spring to ground [m], [N/m]
 * @returns {{x, u, meanDisplacement:number[], walkPerCycle:number[], walkRate:number, endDisplacementStart:number[],
 *   endDisplacementEnd:number[], hotForce:Float64Array, coldForce:Float64Array, anchorForceMin:number,
 *   anchorForceMax:number, iterations:number}} displacements [m], element forces [N] (tension positive)
 */
export function walkingModel({ L, EA, W, mu, steps, cycles = 3, slopeRad = 0, tensionAtEnd = 0, tensionAtStart = 0, uMob = 0.002, anchor, tol = 1e-8 }) {
  const n = steps[0].length, nn = n + 1, h = L / n, kb = EA / h;
  const fLim = mu * W * Math.cos(slopeRad), grav = W * Math.sin(slopeRad);
  const trib = new Float64Array(nn);
  for (let i = 0; i < nn; i++) trib[i] = i === 0 || i === n ? h / 2 : h;
  const u = new Float64Array(nn), up = new Float64Array(nn), q = new Float64Array(nn), res = new Float64Array(nn), slid = new Uint8Array(nn);
  const lower = new Float64Array(nn).fill(-kb), upper = new Float64Array(nn).fill(-kb), diag = new Float64Array(nn);
  const eps = new Float64Array(n), u0 = new Float64Array(nn);
  const iAnchor = anchor ? Math.max(0, Math.min(n, Math.round(anchor.x / h))) : -1;
  const kAnchor = anchor ? anchor.stiffness : 0;
  const scale = fLim * L;
  let iterations = 0;

  // residual r = −∂E/∂u at the current u; returns its max norm
  const residual = () => {
    let norm = 0;
    for (let i = 0; i < nn; i++) {
      const lim = fLim * trib[i], ks = lim / uMob;
      let qi = ks * (u[i] - up[i]);
      if (qi > lim) { qi = lim; slid[i] = 1; } else if (qi < -lim) { qi = -lim; slid[i] = 1; } else slid[i] = 0;
      q[i] = qi;
      let r = grav * trib[i] - qi;
      if (i < n) r += kb * (u[i + 1] - u[i]) - EA * eps[i];
      if (i > 0) r -= kb * (u[i] - u[i - 1]) - EA * eps[i - 1];
      if (i === n) r += tensionAtEnd;
      if (i === 0) r -= tensionAtStart;
      if (i === iAnchor) r -= kAnchor * u[i];
      res[i] = r;
      const a = Math.abs(r); if (a > norm) norm = a;
    }
    return norm;
  };
  const slope = (d) => { let g = 0; for (let i = 0; i < nn; i++) g -= res[i] * d[i]; return g; }; // dE/ds along d

  const equilibrate = () => {
    const maxIter = 300 + 6 * nn;
    for (let it = 0; it < maxIter; it++) {
      const norm = residual();
      if (norm <= tol * scale) break;
      if (it === maxIter - 1) {
        if (norm <= 1e-5 * scale) break; // equilibrium to 0.001 % of the total friction: accept
        throw new Error('The walking model did not reach equilibrium in a load step. Check that friction can hold the pipe against the slope and end tensions, or increase the mobilisation displacement.');
      }
      iterations++;
      for (let i = 0; i < nn; i++) {
        const ks = fLim * trib[i] / uMob;
        diag[i] = (i < n ? kb : 0) + (i > 0 ? kb : 0) + (i === iAnchor ? kAnchor : 0) + (slid[i] ? ks * 1e-9 : ks);
      }
      const d = solveTridiagonal(lower, diag, upper, res);
      u0.set(u);
      // exact line search: the directional derivative is non-decreasing in s, negative at s = 0
      for (let i = 0; i < nn; i++) u[i] = u0[i] + d[i];
      residual();
      if (slope(d) > 0) {
        let lo = 0, hi = 1;
        for (let k = 0; k < 40; k++) {
          const s = 0.5 * (lo + hi);
          for (let i = 0; i < nn; i++) u[i] = u0[i] + s * d[i];
          residual();
          if (slope(d) > 0) hi = s; else lo = s;
        }
        const s = lo > 0 ? lo : hi;
        for (let i = 0; i < nn; i++) u[i] = u0[i] + s * d[i];
      }
      if (!(Math.abs(u[0]) < 1e3 * L) || !(Math.abs(u[n]) < 1e3 * L)) throw new Error('The pipe slides without limit: the slope or end tension exceeds the available axial friction.');
    }
    for (let i = 0; i < nn; i++) up[i] = u[i] - q[i] * uMob / (fLim * trib[i]);
  };
  const mean = () => { let s = 0; for (let i = 0; i < nn; i++) s += u[i] * trib[i]; return s / L; };
  const forces = () => { const f = new Float64Array(n); for (let e = 0; e < n; e++) f[e] = kb * (u[e + 1] - u[e]) - EA * eps[e]; return f; };

  equilibrate(); // as-laid state under the permanent loads
  const meanDisplacement = [], walkPerCycle = [], endStart = [], endEnd = [];
  let prev = mean(), hotForce = forces(), aMin = 0, aMax = 0;
  for (let c = 0; c < cycles; c++) {
    let minSum = Infinity;
    const last = c === cycles - 1;
    if (last && iAnchor >= 0) { aMin = Infinity; aMax = -Infinity; }
    for (const st of steps) {
      eps.set(st);
      equilibrate();
      if (last) {
        const f = forces(); let sum = 0;
        for (let e = 0; e < n; e++) sum += f[e];
        if (sum < minSum) { minSum = sum; hotForce = f; }
        if (iAnchor >= 0) { const fa = kAnchor * u[iAnchor]; aMin = Math.min(aMin, fa); aMax = Math.max(aMax, fa); }
      }
    }
    const m = mean();
    meanDisplacement.push(m); walkPerCycle.push(m - prev); prev = m;
    endStart.push(u[0]); endEnd.push(u[n]);
  }
  const x = Float64Array.from({ length: nn }, (_, i) => i * h);
  return { x, u: Float64Array.from(u), meanDisplacement, walkPerCycle, walkRate: walkPerCycle[walkPerCycle.length - 1], endDisplacementStart: endStart, endDisplacementEnd: endEnd,
    hotForce, coldForce: forces(), anchorForceMin: aMin, anchorForceMax: aMax, iterations };
}
