import test from 'node:test';
import assert from 'node:assert/strict';
import {
  helixGeometry, rectWireSection, roundWireSection, wireStrain, helixLayerStiffness, helixContactPressure,
  wireCurvatureChange, wireBendingStress, helixBendingStiffness, criticalCurvature, wireFrictionStress,
  wireStressAtPosition, momentCurvatureLoop,
} from '../../src/engine/xsection/helix.js';

const close = (a, b, rel = 1e-9, abs = 0) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b) + abs, `${a} vs ${b}`);
const DEG = Math.PI / 180;

test('helix geometry: pitch, wire length and fill fraction by hand', () => {
  const R = 0.1, alpha = 30 * DEG, b = 0.01;
  const g = helixGeometry({ R, alpha, wireWidth: b, nWires: 40 });
  close(g.pitch, 2 * Math.PI * 0.1 / Math.tan(alpha));
  close(g.wireLengthPerPitch, Math.hypot(g.pitch, 2 * Math.PI * R));
  close(g.fillFraction, 40 * 0.01 / (2 * Math.PI * 0.1 * Math.cos(alpha)));
  close(g.curvature, 0.25 / 0.1);
  const g2 = helixGeometry({ R, alpha, wireWidth: b, fillFraction: 0.9 });
  assert.equal(g2.nWires, Math.floor(0.9 * 2 * Math.PI * R * Math.cos(alpha) / b));
  assert.throws(() => helixGeometry({ R, alpha, wireWidth: b, nWires: 80 }), RangeError);
});

test('wire sections: rectangle and circle', () => {
  const r = rectWireSection({ width: 0.012, thickness: 0.004 });
  close(r.A, 48e-6); close(r.In, 0.012 * 0.004 ** 3 / 12); close(r.It, 0.004 * 0.012 ** 3 / 12);
  // square bar: Saint-Venant J = 0.1406 a^4 (Timoshenko & Goodier, Theory of Elasticity, §109)
  close(rectWireSection({ width: 1, thickness: 1 }).J, 0.1406, 2e-3);
  const c = roundWireSection({ diameter: 0.005 });
  close(c.J, Math.PI * 0.005 ** 4 / 32); close(c.In + c.It, c.J);
});

test('helical layer stiffness equals second derivatives of the strain energy from exact helix geometry', () => {
  const n = 46, E = 2.07e11, A = 3.6e-5, R = 0.11, alpha = 33 * DEG;
  for (const lay of [1, -1]) {
    const L = 2 * Math.PI * R / Math.tan(alpha); // one pitch
    const l0 = Math.hypot(L, 2 * Math.PI * R);
    // exact wire length after stretch e, twist t (rad/m) and radial strain x = dR/R
    const U = (e, t, x) => {
      const l = Math.hypot(L * (1 + e), R * (1 + x) * (2 * Math.PI + lay * t * L));
      const ew = l / l0 - 1;
      return n * 0.5 * E * A * ew * ew * l0 / L; // energy per metre of pipe
    };
    const h = [1e-5, 1e-5 / R, 1e-5];
    const d2 = (i, j) => {
      const q = (si, sj) => { const v = [0, 0, 0]; v[i] += si * h[i]; v[j] += sj * h[j]; return U(...v); };
      return (q(1, 1) - q(1, -1) - q(-1, 1) + q(-1, -1)) / (4 * h[i] * h[j]);
    };
    const k = helixLayerStiffness({ n, E, A, R, alpha, lay });
    close(k.kAxial, d2(0, 0), 1e-4);
    close(k.kTorsion, d2(1, 1), 1e-4);
    close(k.kCoupling, d2(0, 1), 1e-4);
    close(k.kAxialRadial, d2(0, 2), 1e-4);
    close(k.kTorsionRadial, d2(1, 2), 1e-4);
    close(k.kRadial, d2(2, 2), 1e-4);
    // closed forms quoted in the literature
    close(k.kAxial, n * E * A * Math.cos(alpha) ** 3);
    close(k.kTorsion, n * E * A * R * R * Math.sin(alpha) ** 2 * Math.cos(alpha));
    close(wireStrain({ alpha, R, lay, epsAxial: 1e-3, twist: 2e-3, dR: 1e-5 }).strain,
      Math.cos(alpha) ** 2 * 1e-3 + lay * R * Math.sin(alpha) * Math.cos(alpha) * 2e-3 + Math.sin(alpha) ** 2 * 1e-5 / R);
  }
});

test('cross-wound symmetric armour pair is torque balanced (zero axial-torsion coupling)', () => {
  const p = { n: 50, E: 2.07e11, A: 4e-5, R: 0.12, alpha: 35 * DEG };
  const a = helixLayerStiffness({ ...p, lay: 1 }), b = helixLayerStiffness({ ...p, lay: -1 });
  assert.equal(a.kCoupling + b.kCoupling, 0);
  assert.ok(a.kCoupling > 0);
  close(a.kAxial + b.kAxial, 2 * p.n * p.E * p.A * Math.cos(p.alpha) ** 3);
});

test('contact pressure equals hoop equilibrium of the wire tensions', () => {
  const n = 40, Fw = 12e3, R = 0.1, alpha = 35 * DEG;
  const q = helixContactPressure({ n, wireForce: Fw, R, alpha });
  // wires crossing a longitudinal cut per metre: n tan(a)/(2 pi R), each with hoop component Fw sin(a)
  const Ntheta = n * Math.tan(alpha) / (2 * Math.PI * R) * Fw * Math.sin(alpha);
  close(q.pressure, Ntheta / R);
  close(q.lineLoadPerWire, Fw * Math.sin(alpha) ** 2 / R);
});

test('stick bending stiffness equals the discrete sum over wires of projected plane-section strains', () => {
  const n = 48, E = 2.07e11, A = 4.8e-5, R = 0.13, alpha = 30 * DEG;
  let EI = 0;
  for (let i = 0; i < n; i++) {
    const y = R * Math.cos(2 * Math.PI * i / n + 0.3);
    const strainPerKappa = Math.cos(alpha) ** 2 * y; // wire strain per unit curvature
    EI += E * A * strainPerKappa * Math.cos(alpha) * y; // axial force component x lever arm
  }
  const sec = rectWireSection({ width: 0.012, thickness: 0.004 });
  const b = helixBendingStiffness({ n, E, A, R, alpha, ...sec });
  close(b.stick, EI, 1e-12);
  close(b.stick, n * E * A * R * R * Math.cos(alpha) ** 3 / 2, 1e-12);
  assert.ok(b.slip > 0 && b.slip < 0.01 * b.stick, 'slip stiffness is orders of magnitude below stick stiffness');
  // straight-wire limit of the slip stiffness: alpha -> 0 gives n E In / 2 + n E It / 2 (wires bend with the pipe)
  const s0 = helixBendingStiffness({ n, E, A, R, alpha: 1e-6, ...sec }).slip;
  close(s0, n * E * (sec.In + sec.It) / 2, 1e-6);
});

/** Numerical curvature of a curve on the bent cylinder (torus) by central differences. */
function curveCurvatures(P, normalAt, z, h) {
  const a = P(z - h), b = P(z), c = P(z + h);
  const d1 = [0, 1, 2].map((i) => (c[i] - a[i]) / (2 * h));
  const d2 = [0, 1, 2].map((i) => (c[i] - 2 * b[i] + a[i]) / (h * h));
  const sp = Math.hypot(...d1);
  const t = d1.map((v) => v / sp);
  const d2t = d2[0] * t[0] + d2[1] * t[1] + d2[2] * t[2];
  const k = d2.map((v, i) => (v - d2t * t[i]) / (sp * sp));
  const nrm = normalAt(z);
  const nxt = [nrm[1] * t[2] - nrm[2] * t[1], nrm[2] * t[0] - nrm[0] * t[2], nrm[0] * t[1] - nrm[1] * t[0]];
  return { kn: -(k[0] * nrm[0] + k[1] * nrm[1] + k[2] * nrm[2]), kg: k[0] * nxt[0] + k[1] * nxt[1] + k[2] * nxt[2], speed: sp };
}

test('loxodromic curvature changes match finite-difference geometry of a helix on a torus', () => {
  const R = 0.1, alpha = 35 * DEG, kappa = 0.02, rho = 1 / kappa, h = 2e-4;
  for (const psiAt of [0, 0.6, Math.PI / 2, 2.2]) {
    const psi = (z) => psiAt + Math.tan(alpha) * z / R; // material line of the cylinder
    const bent = (z) => {
      const f = z / rho, r = rho + R * Math.cos(psi(z));
      return [r * Math.sin(f), R * Math.sin(psi(z)), r * Math.cos(f) - rho];
    };
    const nBent = (z) => { const f = z / rho, p = psi(z); return [Math.cos(p) * Math.sin(f), Math.sin(p), Math.cos(p) * Math.cos(f)]; };
    const straight = (z) => [R * Math.cos(psi(z)), R * Math.sin(psi(z)), z];
    const nStraight = (z) => [Math.cos(psi(z)), Math.sin(psi(z)), 0];
    const kb = curveCurvatures(bent, nBent, 0, h), ks = curveCurvatures(straight, nStraight, 0, h);
    const stretch = kb.speed / ks.speed; // per undeformed wire length
    const ref = wireCurvatureChange({ alpha, kappa, psi: psiAt, curve: 'loxodromic' });
    const tol = 0.02 * kappa; // second-order terms are O(R kappa) = 0.2 % of kappa
    close(kb.kn * stretch - ks.kn, ref.dKn, 0, tol);
    close(Math.abs(kb.kg * stretch - ks.kg), Math.abs(ref.dKt), 0, tol);
    close(stretch - 1, Math.cos(alpha) ** 2 * R * kappa * Math.cos(psiAt), 0, 1e-5);
  }
});

test('geodesic curvature change: Clairaut path on the torus has zero geodesic curvature and dKn = 3 cos^2(a) kappa cos(psi)', () => {
  const R = 0.1, alpha = 35 * DEG, kappa = 0.02, rho = 1 / kappa;
  // path z(psi) from Clairaut: (1 + R kappa cos psi) cos a = cos alpha; dz/dpsi = R cos a / (sin a (1 + R kappa cos psi))
  const dz = (p) => {
    const q = 1 + R * kappa * Math.cos(p), ca = Math.cos(alpha) / q, sa = Math.sqrt(1 - ca * ca);
    return R * ca / (sa * q);
  };
  const zOf = (p) => { // Simpson integration from 0
    const m = 400, hh = p / m; let s = dz(0) + dz(p);
    for (let i = 1; i < m; i++) s += (i % 2 ? 4 : 2) * dz(i * hh);
    return s * hh / 3;
  };
  const P = (p) => { const f = zOf(p) / rho, r = rho + R * Math.cos(p); return [r * Math.sin(f), R * Math.sin(p), r * Math.cos(f) - rho]; };
  const nrm = (p) => { const f = zOf(p) / rho; return [Math.cos(p) * Math.sin(f), Math.sin(p), Math.cos(p) * Math.cos(f)]; };
  for (const psiAt of [0.4, 1.3, 2.5]) {
    const k = curveCurvatures(P, nrm, psiAt, 1e-3);
    const ref = wireCurvatureChange({ alpha, kappa, psi: psiAt, curve: 'geodesic' });
    assert.ok(Math.abs(k.kg) < 1e-4 * kappa + 1e-5, `geodesic curvature ${k.kg}`);
    close(k.kn - Math.sin(alpha) ** 2 / R, ref.dKn, 0, 0.03 * kappa);
    assert.equal(ref.dKt, 0);
  }
});

test('wire bending stress: E t/2 dKn and E b/2 dKt', () => {
  const E = 2.07e11, t = 0.004, b = 0.012, alpha = 35 * DEG, kappa = 0.1;
  const s0 = wireBendingStress({ E, thickness: t, width: b, alpha, kappa, psi: 0 });
  close(s0.sigmaNormal, E * 0.002 * Math.cos(alpha) ** 4 * 0.1);
  close(s0.sigmaTransverse, 0, 0, 1);
  const s90 = wireBendingStress({ E, thickness: t, width: b, alpha, kappa, psi: Math.PI / 2 });
  close(Math.abs(s90.sigmaTransverse), E * 0.006 * (1 + Math.sin(alpha) ** 2) * Math.cos(alpha) * 0.1);
  assert.throws(() => wireCurvatureChange({ alpha, kappa, curve: 'x' }), RangeError);
});

test('critical curvature by hand: f / (E A cos^2 sin), full slip at pi/2 times that', () => {
  const E = 2.0e11, A = 5e-5, R = 0.1, alpha = 30 * DEG, f = 500; // N/m
  const c = criticalCurvature({ E, A, R, alpha, n: 40, frictionPerLength: f });
  // E A = 1.0e7 N; cos^2 30 = 0.75; sin 30 = 0.5 -> kappa_cr = 500 / (1e7 * 0.375) = 1.3333e-4 1/m
  close(c.kappaCr, 1.3333333e-4, 1e-6);
  close(c.kappaFullSlip, Math.PI / 2 * 1.3333333e-4, 1e-6);
  // sigma_f = pi R f /(2 sin a A) = pi * 0.1 * 500 / (2 * 0.5 * 5e-5) = 3.1416 MPa
  close(c.sigmaFrictionMax, Math.PI * 1e6, 1e-9);
  close(c.frictionMoment, 4 / Math.PI * c.momentAtOnset);
  close(c.frictionMoment, 2 / Math.PI * 40 * f * R * R / Math.tan(alpha));
  // pressure form: f = mu (qo + qi) 2 pi R cos(a)/n
  const c2 = criticalCurvature({ E, A, R, alpha, n: 40, mu: 0.1, qOuter: 2e6, qInner: 3e6 });
  close(c2.frictionPerLength, 0.1 * 5e6 * 2 * Math.PI * R * Math.cos(alpha) / 40);
});

test('friction stress: stick below kappa_cr, continuous at the slip boundary, friction moment by quadrature', () => {
  const E = 2.0e11, A = 5e-5, R = 0.1, alpha = 30 * DEG, f = 500, n = 40;
  const c = criticalCurvature({ E, A, R, alpha, n, frictionPerLength: f });
  const a = wireFrictionStress({ E, A, R, alpha, kappa: 0.5 * c.kappaCr, psi: 0.3, frictionPerLength: f });
  assert.equal(a.state, 'stick');
  close(a.sigma, E * 0.75 * R * 0.5 * c.kappaCr * Math.cos(0.3));
  const kap = 1.3 * c.kappaCr;
  const p = wireFrictionStress({ E, A, R, alpha, kappa: kap, psi: 0, frictionPerLength: f });
  assert.equal(p.state, 'partial-slip');
  close((Math.PI / 2 - p.psi0) / Math.cos(p.psi0), 1.3, 1e-9);
  const lo = wireFrictionStress({ E, A, R, alpha, kappa: kap, psi: p.psi0 - 1e-7, frictionPerLength: f });
  const hi = wireFrictionStress({ E, A, R, alpha, kappa: kap, psi: p.psi0 + 1e-7, frictionPerLength: f });
  close(lo.sigma, hi.sigma, 1e-5);
  assert.ok(!lo.slipped && hi.slipped);
  // full slip: moment of the triangular force distribution integrated numerically over the circumference
  let M = 0; const m = 20000;
  for (let i = 0; i < m; i++) {
    const psi = 2 * Math.PI * (i + 0.5) / m;
    const s = wireFrictionStress({ E, A, R, alpha, kappa: 10 * c.kappaCr, psi, frictionPerLength: f }).sigma;
    M += n * s * A * Math.cos(alpha) * R * Math.cos(psi) / m;
  }
  close(M, c.frictionMoment, 1e-6);
  // intrados is the mirror image of the extrados
  close(wireFrictionStress({ E, A, R, alpha, kappa: kap, psi: Math.PI, frictionPerLength: f }).sigma, -p.sigma);
});

test('wire stress recovery superposes mean, friction and local bending', () => {
  const o = { E: 2.07e11, A: 4.8e-5, R: 0.1, alpha: 35 * DEG, thickness: 0.004, width: 0.012, frictionPerLength: 800 };
  const r = wireStressAtPosition({ ...o, kappa: 0.05, psi: 0, sigmaAxial: 300e6 });
  assert.equal(r.state, 'full-slip');
  close(r.sigmaMaxCorner, 300e6 + r.sigmaFriction + Math.abs(r.sigmaNormal) + Math.abs(r.sigmaTransverse));
  close(r.sigmaMaxCorner + r.sigmaMinCorner, 2 * (300e6 + r.sigmaFriction));
});

test('stick -> slip ordering and hysteresis loop area equals 4 Mf (kappa_a - kappa_y)', () => {
  const EIstick = 2.0e6, EIslip = 1.5e4, Mf = 400, ka = 0.02;
  const lp = momentCurvatureLoop({ EIstick, EIslip, frictionMoment: Mf, kappaAmp: ka, nPoints: 2001 });
  assert.ok(EIstick > lp.secantStiffness && lp.secantStiffness > EIslip);
  close(lp.kappaY, Mf / (EIstick - EIslip));
  close(lp.Mmax, Mf + EIslip * ka);
  // enclosed area by the shoelace formula on the returned loop
  let area = 0;
  for (let i = 0; i < lp.kappa.length; i++) {
    const j = (i + 1) % lp.kappa.length;
    area += lp.kappa[i] * lp.moment[j] - lp.kappa[j] * lp.moment[i];
  }
  close(Math.abs(area) / 2, 4 * Mf * (ka - lp.kappaY), 2e-3);
  close(lp.energyPerCycle, 4 * Mf * (ka - lp.kappaY));
  close(lp.dampingRatio, lp.energyPerCycle / (2 * Math.PI * lp.Mmax * ka));
  // below the knee: no dissipation, stick stiffness
  const el = momentCurvatureLoop({ EIstick, EIslip, frictionMoment: Mf, kappaAmp: 0.5 * lp.kappaY });
  assert.equal(el.energyPerCycle, 0);
  close(el.secantStiffness, EIstick);
});

test('invalid input throws RangeError naming the argument', () => {
  assert.throws(() => helixLayerStiffness({ n: 10, E: 2e11, A: 1e-5, R: 0.1, alpha: 2 }), /alpha/);
  assert.throws(() => helixLayerStiffness({ n: 10, E: 2e11, A: 1e-5, R: 0.1, alpha: 0.5, lay: 2 }), /lay/);
  assert.throws(() => criticalCurvature({ E: 2e11, A: -1, R: 0.1, alpha: 0.5, frictionPerLength: 1 }), /A must/);
  assert.throws(() => momentCurvatureLoop({ EIstick: 1, EIslip: 2, frictionMoment: 1, kappaAmp: 1 }), RangeError);
});
