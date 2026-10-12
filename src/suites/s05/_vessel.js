// Vessel-induced hydraulics at a berth: propeller and thruster jets, the return current and water-level
// depression (primary wave) of a passing ship, and ship-generated wash waves. SI units throughout.
//
// Each relation below is either derived here from conservation laws or is an empirical jet / wash relation whose
// coefficients are arguments. The defaults are the values commonly quoted in berth-scour guidance (PIANC WG180,
// 2015; BAW 2010 "Principles for the Design of Bank and Bottom Protection for Inland Waterways") and must be
// confirmed against the edition used on the project.
const G = 9.80665;

/**
 * Propeller (or thruster) jet.
 *
 * Efflux velocity. An ideal actuator at zero advance delivers the power P as kinetic energy flux of a contracted
 * jet of diameter D0:  η P = ½ ρ (π D0²/4) V0³, so V0 = [8 η/(π (D0/Dp)²)]^(1/3) · (P/(ρ Dp²))^(1/3). With the
 * momentum-theory contraction D0 = Dp/√2 and an overall efficiency η of about 0.3 the bracket is 1.15, the value
 * commonly quoted for open propellers; for a ducted propeller or thruster the jet does not contract (D0 = Dp).
 * The coefficient cV0 is an argument: V0 = cV0 (P/(ρ Dp²))^(1/3).
 *
 * Jet decay. Beyond the zone of flow establishment a round turbulent jet spreads linearly and keeps its momentum,
 * which gives V_axis(x) = A·V0·D0/x and a Gaussian profile V(x, r) = V_axis·exp(−b r²/x²). A = 2.8 and b = 15.4
 * are the propeller-jet values of Blaauw & van de Kaa (1978) (arguments cAxis, bRad).
 *
 * Bed velocity. On a bed at distance hp below the axis the velocity is V(x, hp); it is largest where
 * d/dx [x⁻¹ exp(−b hp²/x²)] = 0, i.e. x = √(2b)·hp, giving
 *     Vb,max = A·e^(−½)/√(2b) · V0·D0/hp   (= 0.306 V0 D0/hp for A = 2.8, b = 15.4).
 *
 * P [W] power applied, Dp [m] propeller diameter, hp [m] height of the axis above the bed, ducted [bool].
 * Returns V0, D0, xMax, VbMax [m/s, m], thrust T = ρ (π D0²/4) V0² [N] (momentum flux of the jet) and the
 * functions axis(x) and at(x, r).
 */
export function propellerJet({ P, Dp, hp, rho = 1025, ducted = false, cV0 = 1.15, cAxis = 2.8, bRad = 15.4 }) {
  const D0 = ducted ? Dp : Dp / Math.SQRT2;
  const V0 = cV0 * Math.cbrt(P / (rho * Dp * Dp));
  const axis = (x) => Math.min(V0, (cAxis * V0 * D0) / Math.max(x, 1e-9));
  const at = (x, rr) => axis(x) * Math.exp((-bRad * rr * rr) / Math.max(x * x, 1e-12));
  const xMax = Math.sqrt(2 * bRad) * hp;
  const VbMax = Math.min(V0, ((cAxis * Math.exp(-0.5)) / Math.sqrt(2 * bRad)) * V0 * D0 / hp);
  return { V0, D0, xMax, VbMax, thrust: rho * (Math.PI * D0 * D0 / 4) * V0 * V0, axis, at };
}

/**
 * Return current and water-level depression beside a ship moving at speed V through a waterway of limited
 * cross-section (one-dimensional energy method of Schijf, 1949, without the empirical correction factor).
 * In the frame moving with the ship the flow is steady. With Ac the waterway cross-section, Am the midship
 * section of the ship, bw the water-surface width, Ur the return current and z the depression:
 *     continuity   V·Ac = (V + Ur)·(Ac − Am − bw·z)
 *     Bernoulli    z = [(V + Ur)² − V²]/(2g)
 * Substituting gives one equation in z, solved by bisection on the subcritical branch. Above the limiting speed no
 * subcritical solution exists (the ship cannot exceed it in displacement mode): `limited` is then true and the
 * values at the limit are returned. Returns { Ur [m/s], z [m], limited, Vlim [m/s] }.
 */
export function returnCurrent({ V, Ac, Am, bw }) {
  if (!(Ac > Am)) throw new Error('The waterway cross-section must exceed the midship section of the passing vessel.');
  // residual f(z) = V·Ac − √(V² + 2 g z)·(Ac − Am − bw z); f(0) = V·Am > 0, and a root exists only below the limit speed
  const f = (v, z) => v * Ac - Math.sqrt(v * v + 2 * G * z) * (Ac - Am - bw * z);
  const zTop = (Ac - Am) / bw;
  const solve = (v) => {
    // locate the first sign change from z = 0
    let zPrev = 0; const n = 400;
    for (let i = 1; i <= n; i++) {
      const z = (zTop * i) / n;
      if (f(v, z) <= 0) { let a = zPrev, b = z; for (let k = 0; k < 60; k++) { const m = 0.5 * (a + b); if (f(v, m) > 0) a = m; else b = m; } return 0.5 * (a + b); }
      zPrev = z;
    }
    return null;
  };
  let z = solve(V), limited = false, Vlim = NaN;
  // limiting speed: the largest V with a solution
  { let a = 0.01, b = Math.sqrt(G * Ac / bw); for (let k = 0; k < 40; k++) { const m = 0.5 * (a + b); if (solve(m) !== null) a = m; else b = m; } Vlim = a; }
  if (z === null) { limited = true; z = solve(Vlim * 0.999) ?? 0; V = Vlim * 0.999; }
  return { Ur: Math.sqrt(V * V + 2 * G * z) - V, z, limited, Vlim };
}

/**
 * Height of the secondary (wash) waves of a moving vessel at distance s from its track, in water of depth h:
 *     H = α·h·(s/h)^(−1/3)·Fr_h⁴,   Fr_h = V/√(g h)
 * The cube-root decay with distance follows from the dispersion of the diverging Kelvin wave system; the
 * fourth-power dependence on the depth Froude number and the coefficient α are empirical (α about 1.0 for tugs,
 * patrol boats and loaded conventional vessels, 0.5 for empty push units, 0.35 for empty conventional vessels in
 * the Delft measurements reported in PIANC and the Rock Manual). α is an argument. Valid for Fr_h below about 0.8.
 * The period of the largest waves follows from the Kelvin angle: the wave travels at V·cos 35.3°, so in deep
 * water T = 2π·0.816·V/g = 5.1 V/g.
 */
export function washWave({ V, h, s, alpha = 1 }) {
  const Fr = V / Math.sqrt(G * h);
  return { H: alpha * h * Math.cbrt(h / Math.max(s, 0.1 * h)) * Fr ** 4, T: (2 * Math.PI * Math.cos((35.27 * Math.PI) / 180) * V) / G, Fr };
}
