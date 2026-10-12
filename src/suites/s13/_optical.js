// Fibre helpers of the optical-link module that replace an engine function with a defect.
// The engine's refractiveIndexProfile rejects its own documented default exponent (g = Infinity, the step profile)
// because the argument check requires a finite number; the power-law profile is therefore evaluated here.

/**
 * Power-law (α-profile) refractive index n(r) = n1·sqrt(1 − 2Δ·(r/a)^g) inside the core and n1·sqrt(1 − 2Δ) in the
 * cladding; g = Infinity is the step profile, g = 2 the parabolic profile.
 * Reference: Gloge & Marcatili, Bell Syst. Tech. J. 52 (1973) 1563.
 * @param {object} o r [m], coreRadius [m], n1 index on the axis, delta relative index difference, g profile exponent
 * @returns {number} refractive index at r
 */
export function indexProfile({ r, coreRadius, n1, delta, g = Infinity }) {
  if (!(r >= 0) || !(coreRadius > 0) || !(n1 > 0)) throw new Error('Refractive-index profile: radius, core radius and core index must be positive.');
  if (!(delta > 0 && delta < 0.5)) throw new Error('Refractive-index profile: the relative index difference must lie between 0 and 0.5.');
  if (!(g > 0)) throw new Error('Refractive-index profile: the profile exponent must be positive.');
  if (r >= coreRadius) return n1 * Math.sqrt(1 - 2 * delta);
  const x = g === Infinity ? 0 : (r / coreRadius) ** g;
  return n1 * Math.sqrt(1 - 2 * delta * x);
}
