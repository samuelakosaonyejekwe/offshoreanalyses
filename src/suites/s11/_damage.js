// Damaged floating equilibrium by physical relaxation, used when the Newton solver of the hydrostatics engine
// stalls on non-smooth cases (a flooded compartment whose top or bottom crosses the waterline makes the
// restoring moment piecewise). Heel and trim follow the residual levers between B and G with step halving on
// every sign change, the draft is corrected with the net waterplane area: this always walks to a stable
// equilibrium if one exists within the capsize bound.
import { bodyIntegrals, bodyToEarth, boxMesh } from '../../engine/hydro/hydrostatics.js';

export const DEG = Math.PI / 180;

/** Compartment records { name, box, permeability } → lost-buoyancy list with meshes. */
export function lostFrom(rows) {
  return rows.map((r) => ({ name: r.name, mesh: boxMesh(r.box), permeability: r.permeability, box: r.box }));
}

/**
 * @param {{mesh, lost, mass (kg), cog [x,y,z] body, rho, gmT, gmL (m, positive scales), draft0, maxDeg}} o
 * @returns {{draft, heel, trim (rad), converged, capsized, iterations, leverX, leverY (m)}}
 * Frame convention of the engine: p_earth = Ry(trim)·Rx(heel)·p_body − (0, 0, draft); buoyancy acting at an earth-frame
 * offset (dx, dy) from G gives moments +dy about x (increasing heel) and −dx about y (decreasing trim).
 */
export function damagedEquilibrium({ mesh, lost, mass, cog, rho = 1025, gmT, gmL, draft0, maxDeg = 60 }) {
  const V0 = mass / rho;
  let T = draft0; let h = 0; let t = 0;
  let ah = 0.6; let at = 0.6; let sh = 0; let st = 0;
  const sT = Math.max(Math.abs(gmT), 0.5); const sL = Math.max(Math.abs(gmL), 0.5);
  const cap = 1.5 * DEG;
  let dy = 0; let dx = 0; let it = 0; let converged = false; let capsized = false;
  for (; it < 600; it++) {
    let r;
    // draft: Newton on volume with the net waterplane area (3 passes are ample between small attitude steps)
    for (let k = 0; k < 4; k++) {
      r = bodyIntegrals(mesh, lost, T, h, t);
      const Aw = Math.max(r.Aw, 1e-3 * V0 ** (2 / 3));
      const dT = (V0 - r.V) / Aw;
      T += Math.max(-2, Math.min(2, dT));
      if (Math.abs(dT) < 1e-7) break;
    }
    r = bodyIntegrals(mesh, lost, T, h, t);
    if (!(r.V > 0)) { capsized = true; break; }
    const g = bodyToEarth(cog, T, h, t);
    dy = r.My / r.V - g[1]; dx = r.Mx / r.V - g[0];
    if (Math.abs(dy) < 2e-6 * sT && Math.abs(dx) < 2e-6 * sL && Math.abs(V0 - r.V) < 1e-7 * V0) { converged = true; break; }
    if (sh && Math.sign(dy) !== sh) ah *= 0.5;
    if (st && Math.sign(dx) !== st) at *= 0.5;
    sh = Math.sign(dy); st = Math.sign(dx);
    h += Math.max(-cap, Math.min(cap, ah * dy / sT));
    t -= Math.max(-cap, Math.min(cap, at * dx / sL));
    if (Math.abs(h) > maxDeg * DEG || Math.abs(t) > maxDeg * DEG) { capsized = true; break; }
  }
  return { draft: T, heel: h, trim: t, converged, capsized, iterations: it, leverX: dx, leverY: dy };
}
