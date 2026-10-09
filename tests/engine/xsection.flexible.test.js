import test from 'node:test';
import assert from 'node:assert/strict';
import {
  buildFlexibleStack, axisymmetricResponse, layerUtilisations, burstPressure, axialFailureTension, torsionLimits,
  bendingProperties, minimumBendRadius, massAndWeight, equivalentBeamProperties, thermalProfile,
  pressureThermalElongation, exampleFlexibleDesign, carcassBucklingFactor, API17J_FACTORS, SHEATH_STRAIN_LIMITS,
} from '../../src/engine/xsection/flexible.js';

const close = (a, b, rel = 1e-9, abs = 0) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b) + abs, `${a} vs ${b}`);
const DEG = Math.PI / 180;
const riser = () => buildFlexibleStack(exampleFlexibleDesign('riser-6in'));

/** Independent check of statics and kinematics of a solution, written from the governing equations. */
function checkSolution(stack, r, { tension = 0, pInt = 0, pExt = 0, pAnn = 0, torque = 0 }) {
  const Ls = stack.layers, s1 = Ls[stack.sealInner], s2 = Ls[stack.sealOuter];
  const fScale = Math.max(Math.abs(tension), pInt * Math.PI * s1.rIn ** 2, pExt * Math.PI * s2.rOut ** 2, 1);
  const pScale = Math.max(pInt, pExt, fScale / (Math.PI * s2.rOut ** 2));
  // global axial equilibrium: wall force = effective tension + end caps on the sealing layers
  const endCap = s1 === s2 ? Math.PI * (pInt * s1.rIn ** 2 - pExt * s1.rOut ** 2)
    : Math.PI * (pInt * s1.rIn ** 2 - pAnn * s1.rOut ** 2 + pAnn * s2.rIn ** 2 - pExt * s2.rOut ** 2);
  let sumN = 0, sumM = 0, sumHoop = 0;
  Ls.forEach((L, i) => {
    const q = r.layers[i];
    sumN += q.axialForce; sumM += q.torque; sumHoop += q.hoopForce;
    // thin-ring statics of each layer
    close(q.hoopForce, q.pIn * L.rIn - q.pOut * L.rOut, 1e-7, 1e-9 * pScale * L.R);
    assert.ok(q.contactIn >= -1e-9 * pScale && q.contactOut >= -1e-9 * pScale, 'contact pressures are compressive');
    if (L.kind === 'helical') {
      const c = Math.cos(L.alpha), s = Math.sin(L.alpha);
      const ew = c * c * r.axialStrain + L.lay * L.R * s * c * r.twist + s * s * q.radialDisp / L.R;
      close(q.wireStress, L.E * ew, 1e-7, 1e-9 * pScale);
      close(q.axialForce, q.wireStress * L.nA * c, 1e-9, 1e-9 * fScale);
      close(q.torque, L.lay * q.wireStress * L.nA * L.R * s, 1e-9, 1e-9 * fScale);
      close(q.hoopForce, q.wireStress * L.nA * s * s / (2 * Math.PI * L.R * c), 1e-9, 1e-9 * pScale * L.R);
    } else {
      // 3-D Hooke's law inverted: strains from the reported stresses
      const ez = (q.sigmaAxial - L.nu * (q.sigmaHoop + q.sigmaRadial)) / L.E;
      const et = (q.sigmaHoop - L.nu * (q.sigmaAxial + q.sigmaRadial)) / L.E;
      const er = (q.sigmaRadial - L.nu * (q.sigmaAxial + q.sigmaHoop)) / L.E;
      const eScale = pScale / L.E;
      close(ez, r.axialStrain, 1e-7, 1e-9 * eScale);
      close(et, q.radialDisp / L.R, 1e-7, 1e-9 * eScale);
      close(er, q.thicknessChange / L.t, 1e-7, 1e-9 * eScale);
      close(q.sigmaRadial * L.R, -(q.pIn * L.rIn + q.pOut * L.rOut) / 2, 1e-7, 1e-9 * pScale * L.R);
      close(q.axialForce, q.sigmaAxial * 2 * Math.PI * L.R * L.t, 1e-9, 1e-9 * fScale);
    }
  });
  close(sumN, tension + endCap, 1e-8, 1e-9 * fScale);
  close(sumM, torque, 1e-8, 1e-9 * fScale * s2.rOut);
  // global ring statics: contact pressures cancel, only fluid pressures on the seals remain
  const hoopExpected = s1 === s2 ? pInt * s1.rIn - pExt * s1.rOut : pInt * s1.rIn - pAnn * s1.rOut + pAnn * s2.rIn - pExt * s2.rOut;
  close(sumHoop, hoopExpected, 1e-7, 1e-9 * pScale * s2.rOut);
  // compatibility at closed interfaces, non-negative gap at open ones
  r.interfaces.forEach((f, j) => {
    const a = r.layers[j], b = r.layers[j + 1];
    const gap = (b.radialDisp - b.thicknessChange / 2) - (a.radialDisp + a.thicknessChange / 2);
    if (f.closed) assert.ok(Math.abs(gap) < 1e-12 * s2.rOut + 1e-9 * Math.abs(a.radialDisp), `compatibility residual ${gap}`);
    else { assert.ok(gap >= 0); assert.equal(f.contactPressure, 0); }
  });
}

test('stack geometry: radii accumulate, metal area from fill fraction', () => {
  const st = riser();
  let r = 0.0762;
  for (const L of st.layers) { close(L.rIn, r, 1e-12); r += L.t; close(L.R, L.rIn + L.t / 2, 1e-12); }
  close(st.outerRadius, r, 1e-12);
  const ta = st.layers[4];
  close(ta.nA, 0.9 * 2 * Math.PI * ta.R * ta.t * Math.cos(35 * DEG), 1e-12);
  assert.equal(st.layers[st.sealInner].type, 'pressureSheath');
  assert.equal(st.layers[st.sealOuter].type, 'outerSheath');
  assert.throws(() => buildFlexibleStack({ innerDiameter: 0.1, layers: [{ thickness: 0.01, E: 1e9, alphaDeg: 30, fill: 0.9 }] }), /sealing/);
});

test('single isotropic tube in tension: eps = T/(EA), hoop strain = -nu eps, thickness strain = -nu eps', () => {
  const st = buildFlexibleStack({ innerDiameter: 0.2, layers: [{ type: 'pressureSheath', thickness: 0.01, E: 1e9, nu: 0.4 }] });
  const T = 5e4, A = 2 * Math.PI * 0.105 * 0.01;
  const r = axisymmetricResponse(st, { tension: T });
  close(r.axialStrain, T / (1e9 * A), 1e-9);
  close(r.layers[0].radialDisp / 0.105, -0.4 * r.axialStrain, 1e-9);
  close(r.layers[0].thicknessChange / 0.01, -0.4 * r.axialStrain, 1e-9);
  close(r.stiffness.K11, 1e9 * A, 1e-9);
  close(r.stiffness.K22, 1e9 / 2.8 * Math.PI * (0.11 ** 4 - 0.1 ** 4) / 2, 1e-9);
});

test('single tube, internal pressure with closed ends: thin-ring hoop force p r_i and end-cap axial force p pi r_i^2', () => {
  const st = buildFlexibleStack({ innerDiameter: 0.2, layers: [{ type: 'pressureSheath', thickness: 0.004, E: 2e11, nu: 0.3 }] });
  const p = 5e6;
  const r = axisymmetricResponse(st, { pInt: p });
  close(r.layers[0].sigmaHoop, p * 0.1 / 0.004, 1e-9);
  close(r.layers[0].sigmaAxial, p * Math.PI * 0.01 / (2 * Math.PI * 0.102 * 0.004), 1e-9);
  close(r.tensionEffective, 0, 0, 1e-6);
  checkSolution(st, r, { pInt: p });
});

test('riser: equilibrium, constitutive and compatibility residuals vanish for combined loads', () => {
  const st = riser();
  const cases = [
    { tension: 5e5 },
    { pInt: 2e7 },
    { tension: 3e5, pInt: 2e7, pExt: 1e7, torque: 2e3 },
    { tension: 8e5, pInt: 1e7, pExt: 1.5e7, torque: -3e3 },
    { pExt: 1e7 },
  ];
  for (const c of cases) {
    const r = axisymmetricResponse(st, c);
    checkSolution(st, r, c);
    assert.ok(Math.abs(r.residuals.axial) < 1e-6 * (Math.abs(r.wallForce) + 1));
    close(r.tensionEffective, c.tension ?? 0, 1e-8, 1e-3);
  }
});

test('riser, pure internal pressure with closed ends: pressure armour carries the hoop load, tensile armour the end cap', () => {
  const st = riser();
  const p = 2e7;
  const r = axisymmetricResponse(st, { pInt: p });
  const ps = st.layers[st.sealInner];
  const hoopTotal = p * ps.rIn; // N per metre, ring statics of the whole wall
  const endCap = p * Math.PI * ps.rIn ** 2;
  const pa = r.layers[2], ta1 = r.layers[4], ta2 = r.layers[6];
  let sumHoop = 0, sumAx = 0;
  for (const q of r.layers) { sumHoop += q.hoopForce; sumAx += q.axialForce; }
  close(sumHoop, hoopTotal, 1e-8);
  close(sumAx, endCap, 1e-8);
  assert.ok(pa.hoopForce > 0.6 * hoopTotal, 'pressure armour takes most of the hoop load');
  assert.ok(ta1.axialForce + ta2.axialForce > 0.9 * endCap, 'tensile armour takes the end-cap load');
  // carcass sits in the bore fluid: unloaded apart from the small contact with the expanding sheath
  assert.ok(Math.abs(r.layers[0].wireStress) < 0.02 * pa.wireStress);
  // pressure armour hoop stress by thin-ring statics from its own face pressures
  const L = st.layers[2];
  close(pa.wireStress * L.nA * Math.sin(L.alpha) ** 2 / (2 * Math.PI * L.R * Math.cos(L.alpha)), pa.pIn * L.rIn - pa.pOut * L.rOut, 1e-8);
});

test('tangent stiffness matrix is symmetric and reproduces the solved strain and twist', () => {
  const st = riser();
  const c = { tension: 6e5, pInt: 1.5e7, torque: 1.5e3 };
  const r = axisymmetricResponse(st, c);
  const K = r.stiffness;
  close(K.K12, K.K21, 1e-7, 1e-9 * Math.sqrt(K.K11 * K.K22));
  // incremental check: an extra tension dT changes (eps, tau) by K^-1 [dT, 0]
  const dT = 1e3;
  const r2 = axisymmetricResponse(st, { ...c, tension: c.tension + dT });
  const det = K.K11 * K.K22 - K.K12 * K.K21;
  close(r2.axialStrain - r.axialStrain, K.K22 * dT / det, 1e-6);
  close(r2.twist - r.twist, -K.K21 * dT / det, 1e-6, 1e-12);
  // cross-wound armour: coupling is weak compared with sqrt(K11 K22)
  assert.ok(Math.abs(K.K12) < 0.15 * Math.sqrt(K.K11 * K.K22));
});

test('two identical cross-wound armour layers at the same lay angle: tension produces almost no twist', () => {
  const st = riser();
  const r = axisymmetricResponse(st, { tension: 1e6 });
  // the two layers differ only in radius (4 mm armour + tape) and hence wire count, so the residual twist is a
  // small fraction of the twist that would unload a single layer (eps cos^2 + R sin cos tau = 0)
  const twistIfOneLayer = r.axialStrain / (st.layers[4].R * Math.tan(35 * DEG));
  assert.ok(Math.abs(r.twist) < 0.1 * twistIfOneLayer);
  assert.ok(r.layers[4].torque * r.layers[6].torque < 0, 'layer torques oppose each other');
});

test('prescribed axial strain and twist reproduce the force-controlled solution', () => {
  const st = riser();
  const c = { tension: 4e5, pInt: 1e7, pExt: 5e6, torque: 1e3 };
  const r = axisymmetricResponse(st, c);
  const r2 = axisymmetricResponse(st, { pInt: 1e7, pExt: 5e6, axialStrain: r.axialStrain, twist: r.twist });
  close(r2.tensionEffective, c.tension, 1e-7);
  close(r2.torque, c.torque, 1e-6);
});

test('active-set contact: lay-opening torque opens an armour interface, pressures stay non-negative', () => {
  const st = riser();
  const c = { tension: 1e4, torque: -4e4 };
  const r = axisymmetricResponse(st, c);
  assert.ok(r.interfaces.some((f) => !f.closed), 'an interface opened');
  checkSolution(st, r, c);
  for (const f of r.interfaces) { assert.ok(f.contactPressure >= 0); assert.ok(f.gap >= 0); }
});

test('flooded annulus: a free outer sheath is in a hydrostatic state, the pressure sheath carries the external pressure', () => {
  // two sealing shells separated by an initial gap: the outer one sees the same pressure on both faces
  const st2 = buildFlexibleStack({ innerDiameter: 0.2, layers: [
    { type: 'pressureSheath', thickness: 0.01, E: 1e9, nu: 0.4 },
    { type: 'outerSheath', thickness: 0.01, E: 1e9, nu: 0.4 },
  ] });
  const p = 1.2e7;
  const h = axisymmetricResponse(st2, { pExt: p, annulus: 'flooded', initialGaps: [5e-3] });
  // thin-shell statics give exactly -p in the hoop and axial directions and -p (r_in + r_out)/(2R) = -p radially
  close(h.layers[1].sigmaHoop, -p, 1e-9);
  close(h.layers[1].sigmaRadial, -p, 1e-9);
  close(h.layers[1].sigmaAxial, h.layers[1].axialForce / (2 * Math.PI * 0.115 * 0.01), 1e-9);
  assert.equal(h.interfaces[0].closed, false);
  const st = riser();
  const r = axisymmetricResponse(st, { pExt: p, annulus: 'flooded' });
  assert.ok(r.layers[0].wireStress < 0, 'carcass in hoop compression');
  assert.ok(r.layers[1].pOut >= p, 'annulus water acts on the pressure sheath');
  checkSolution(st, r, { pExt: p, pAnn: p });
  // dry annulus: external pressure reaches the pressure armour and carcass through the armours
  const dry = axisymmetricResponse(st, { pExt: p });
  assert.ok(dry.layers[2].wireStress < 0);
  checkSolution(st, dry, { pExt: p });
});

test('burst and axial failure: governing layer reaches its strength exactly at the reported load', () => {
  const st = riser();
  const b = burstPressure(st);
  assert.equal(b.burst, Math.min(b.pressureArmour, b.tensileArmour));
  const r = axisymmetricResponse(st, { pInt: b.pressureArmour });
  close(r.layers[2].wireStress, 1000e6, 1e-6);
  const rt = axisymmetricResponse(st, { pInt: b.tensileArmour });
  close(Math.max(rt.layers[4].wireStress, rt.layers[6].wireStress), 1400e6, 1e-6);
  // order of magnitude by thin-ring statics: p ~ sigma_u x (metal hoop thickness) / r
  const L = st.layers[2];
  const pRing = 1000e6 * L.nA * Math.sin(L.alpha) ** 2 / (2 * Math.PI * L.R * Math.cos(L.alpha)) / st.layers[1].rIn;
  assert.ok(b.pressureArmour > pRing && b.pressureArmour < 1.6 * pRing);
  const a = axialFailureTension(st);
  const ra = axisymmetricResponse(st, { tension: a.tension });
  close(Math.max(ra.layers[4].wireStress, ra.layers[6].wireStress), 1400e6, 1e-6);
  // upper bound: all armour metal at UTS along the wire, projected on the axis
  const cap = (st.layers[4].nA + st.layers[6].nA) * 1400e6 * Math.cos(35 * DEG);
  assert.ok(a.tension < 1.05 * cap && a.tension > 0.8 * cap);
});

test('utilisation factors and torsion limits', () => {
  const st = riser();
  const r = axisymmetricResponse(st, { tension: 5e5, pInt: 2e7 });
  const u = layerUtilisations(st, r, { condition: 'normal' });
  const pa = u.layers.find((x) => x.type === 'pressureArmour');
  close(pa.utilisation, Math.abs(r.layers[2].wireStress) / (0.55 * 0.9 * 1000e6));
  const h = layerUtilisations(st, r, { condition: 'hydrotest' });
  close(h.layers.find((x) => x.type === 'tensileArmour').factor, 0.91);
  assert.equal(API17J_FACTORS.tensileArmour.normal, 0.67);
  close(carcassBucklingFactor(600).factor, 0.76);
  assert.equal(carcassBucklingFactor(100).factor, 0.67);
  assert.equal(carcassBucklingFactor(2000).factor, 0.85);
  const t = torsionLimits(st, { loads: { tension: 3e5, pInt: 1e7 } });
  assert.ok(t.tightening.governing * st.layers[6].lay > 0);
  assert.ok(t.opening.governing * st.layers[6].lay < 0);
  // at the stress-limited tightening torque the worst armour wire sits at 0.67 x 0.9 UTS
  const rs = axisymmetricResponse(st, { tension: 3e5, pInt: 1e7, torque: t.tightening.stressLimit });
  close(Math.max(Math.abs(rs.layers[4].wireStress), Math.abs(rs.layers[6].wireStress)), 0.67 * 0.9 * 1400e6, 1e-6);
  assert.ok(Math.abs(t.opening.separation) <= Math.abs(t.opening.stressLimit ?? Infinity));
});

test('bending: stick > slip stiffness, critical curvature of an armour layer by hand from its contact pressures', () => {
  const st = riser();
  const r = axisymmetricResponse(st, { tension: 5e5, pInt: 2e7 });
  const mu = 0.12;
  const b = bendingProperties(st, { response: r, mu });
  assert.ok(b.EIstick > 10 * b.EIslip && b.EIslip > 0);
  const L = st.layers[6], q = r.layers[6];
  // f/A = mu q 2 pi R cos(a)/(n A) = mu q/(F t) with n A = F 2 pi R t cos(a); kappa_cr = (f/A)/(E cos^2(a) sin(a))
  const byHand = mu * (q.contactIn + q.contactOut) / (L.fill * L.t) / (L.E * Math.cos(L.alpha) ** 2 * Math.sin(L.alpha));
  close(b.layers[1].kappaCr, byHand, 1e-9);
  close(b.kappaCrMin, Math.min(b.layers[0].kappaCr, b.layers[1].kappaCr));
  // polymer-only part of the slip stiffness
  let poly = 0;
  for (const s of st.layers) if (s.kind === 'shell') poly += s.E * Math.PI * (s.rOut ** 4 - s.rIn ** 4) / 4;
  assert.ok(b.EIslip > poly);
  const eq = equivalentBeamProperties(st, { loads: { tension: 5e5, pInt: 2e7 }, mu });
  close(eq.EIstick, b.EIstick); close(eq.kappaCr, b.kappaCrMin);
  assert.ok(eq.EA < eq.EAfixed * (1 + 1e-12) && eq.GJ <= eq.GJfixed * (1 + 1e-12));
});

test('minimum bend radius: locking geometry and sheath strain by hand', () => {
  const st = riser();
  const m = minimumBendRadius(st);
  const carc = st.layers[0], pa = st.layers[2];
  const lock = Math.max(carc.rOut / 0.09, pa.rOut / 0.07);
  close(m.lockingRadius, lock);
  const sheath = st.layers[7].rOut / 0.077;
  close(m.sheathStaticRadius, sheath);
  close(m.storage, Math.max(1.1 * lock, sheath));
  close(m.dynamic, 1.5 * m.storage);
  close(minimumBendRadius(st, { storageFactor: 1, dynamicFactor: 2 }).dynamic, 2 * Math.max(lock, sheath));
  assert.equal(SHEATH_STRAIN_LIMITS.PVDF.dynamic, 0.035);
});

test('mass and weight by hand for a two-layer pipe', () => {
  const st = buildFlexibleStack({ innerDiameter: 0.2, layers: [
    { type: 'pressureSheath', thickness: 0.01, E: 1e9, nu: 0.4, rho: 1000 },
    { type: 'tensileArmour', thickness: 0.005, alphaDeg: 30, lay: 1, fill: 0.8, E: 2e11, rho: 7850 },
    { type: 'outerSheath', thickness: 0.01, E: 1e9, nu: 0.4, rho: 950 },
  ] });
  const mSheath = 1000 * Math.PI * (0.11 ** 2 - 0.1 ** 2);
  const mArm = 7850 * 0.8 * Math.PI * (0.115 ** 2 - 0.11 ** 2); // fill x annulus area
  const mOuter = 950 * Math.PI * (0.125 ** 2 - 0.115 ** 2);
  const w = massAndWeight(st, { rhoContent: 800, rhoSeawater: 1025, g: 9.81 });
  close(w.massEmpty, mSheath + mArm + mOuter, 1e-12);
  close(w.contentArea, Math.PI * 0.01, 1e-12);
  close(w.massOperating, w.massEmpty + 800 * Math.PI * 0.01, 1e-12);
  close(w.submergedEmpty, (w.massEmpty - 1025 * Math.PI * 0.125 ** 2) * 9.81, 1e-12);
  const fl = massAndWeight(st, { annulusFlooded: true });
  close(fl.annulusWaterMass, 1025 * 0.2 * Math.PI * (0.115 ** 2 - 0.11 ** 2), 1e-9);
});

test('thermal profile: composite-cylinder resistances and U-value by hand', () => {
  const st = buildFlexibleStack({ innerDiameter: 0.2, layers: [
    { type: 'pressureSheath', thickness: 0.01, E: 1e9, k: 0.3 },
    { type: 'outerSheath', thickness: 0.02, E: 1e9, k: 0.4 },
  ] });
  const R1 = Math.log(0.11 / 0.1) / (2 * Math.PI * 0.3), R2 = Math.log(0.13 / 0.11) / (2 * Math.PI * 0.4);
  const Rin = 1 / (2 * Math.PI * 0.1 * 500), Rout = 1 / (2 * Math.PI * 0.13 * 1000);
  const t = thermalProfile(st, { TinC: 80, ToutC: 4, hIn: 500, hOut: 1000 });
  const q = 76 / (R1 + R2 + Rin + Rout);
  close(t.heatLoss, q);
  close(t.UvalueInner, q / (2 * Math.PI * 0.1 * 76));
  close(t.interfaceTempsC[1], 80 - q * (Rin + R1));
  close(t.interfaceTempsC[2], 4 + q * Rout);
});

test('pressure and thermal elongation: free single tube expands by alphaT dT; closed-end thin tube under pressure', () => {
  const st = buildFlexibleStack({ innerDiameter: 0.2, layers: [{ type: 'pressureSheath', thickness: 0.002, E: 2e11, nu: 0.3, alphaT: 1.2e-5 }] });
  const e = pressureThermalElongation(st, { pInt: 1e6, deltaT: 50 });
  close(e.strainThermal, 1.2e-5 * 50, 1e-9);
  // thin shell with through-thickness stress: eps_z = [sz - nu (st + sr)]/E from the statics of the model
  const R = 0.101, t = 0.002, p = 1e6;
  const sz = p * Math.PI * 0.01 / (2 * Math.PI * R * t), sth = p * 0.1 / t, sr = -p * 0.1 / (2 * R);
  close(e.strainPressure, (sz - 0.3 * (sth + sr)) / 2e11, 1e-9);
  close(e.strainTotal, e.strainPressure + e.strainThermal, 1e-9);
  // the riser lengthens or shortens only slightly under pressure but always expands when heated
  const er = pressureThermalElongation(riser(), { pInt: 2e7, deltaT: 40 });
  assert.ok(er.strainThermal > 0 && Number.isFinite(er.strainPressure));
});

test('example designs are labelled generic and solvable', () => {
  for (const id of ['riser-6in', 'flowline-10in']) {
    const d = exampleFlexibleDesign(id);
    assert.equal(d.generic, true);
    assert.match(d.label, /demonstration/);
    const st = buildFlexibleStack(d);
    const r = axisymmetricResponse(st, { tension: 2e5, pInt: 1e7, pExt: 3e6 });
    checkSolution(st, r, { tension: 2e5, pInt: 1e7, pExt: 3e6 });
  }
  assert.throws(() => exampleFlexibleDesign('x'), RangeError);
  assert.throws(() => axisymmetricResponse(riser(), { pInt: -1 }), /pInt/);
});
