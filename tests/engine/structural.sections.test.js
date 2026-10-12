import test from 'node:test';
import assert from 'node:assert/strict';
import {
  tubeSection, solidRoundSection, boxSection, rectHollowSection, rectSection, iSection, rolledISection, iSectionByName, EURO_I_PROFILES,
  weldedISection, plateGirderSection, channelSection, angleSection, teeSection, stiffenedPlateStrip, polygonSection, compositeTubeSection,
  OFFSHORE_TUBULARS, tubularByName, resolveSection,
} from '../../src/engine/structural/sections.js';

const near = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);
const ring = (r, n) => Array.from({ length: n }, (_, i) => [r * Math.cos(2 * Math.PI * i / n), r * Math.sin(2 * Math.PI * i / n)]);

test('tube: hand values and agreement with the polygon integrator', () => {
  const s = tubeSection({ D: 0.5, t: 0.02 });
  // A = π/4 (0.5² − 0.46²) = 0.0301593; I = π/64 (0.5⁴ − 0.46⁴) = 0.0490874 × 0.01772544 = 8.700955e-4
  near(s.A, 0.0301593, 1e-6);
  near(s.Iy, 8.700955e-4, 1e-6);
  near(s.J, 2 * s.Iy, 1e-15);
  near(s.Wply, (0.5 ** 3 - 0.46 ** 3) / 6, 1e-15);
  near(s.iy, Math.sqrt(0.5 ** 2 + 0.46 ** 2) / 4, 1e-12, 'r = √(D² + d²)/4');
  const p = polygonSection({ outer: ring(0.25, 2000), holes: [ring(0.23, 2000)] });
  near(p.A, s.A, 1e-5); near(p.Iy, s.Iy, 1e-5); near(p.Iz, s.Iz, 1e-5);
  // Cowper shear coefficients: solid circle 6(1+ν)/(7+6ν) = 0.8864; thin tube → 2(1+ν)/(4+3ν) = 0.5306
  near(solidRoundSection({ D: 0.2 }).Asy / (Math.PI * 0.01), 0.886364, 1e-6);
  const thin = tubeSection({ D: 1, t: 1e-5 });
  near(thin.Asy / thin.A, 2 * 1.3 / 4.9, 1e-4);
  assert.throws(() => tubeSection({ D: 0.5, t: 0.3 }), RangeError);
});

test('rectangle, box and RHS', () => {
  const r = rectSection({ h: 0.4, b: 0.2 });
  near(r.Iy, 0.2 * 0.4 ** 3 / 12, 1e-14);
  near(r.Wply, 0.2 * 0.16 / 4, 1e-14);
  // torsion constant of a 2:1 rectangle: 0.229 a b³ (Roark)
  near(r.J, 0.229 * 0.4 * 0.2 ** 3, 2e-3);
  near(rectSection({ h: 0.1, b: 0.1 }).J, 0.1406 * 1e-4, 2e-3, 'square: 0.1406 a⁴');
  const b = boxSection({ h: 0.6, b: 0.4, tw: 0.02, tf: 0.03 });
  const p = polygonSection({ outer: [[-0.2, -0.3], [0.2, -0.3], [0.2, 0.3], [-0.2, 0.3]], holes: [[[-0.18, -0.27], [0.18, -0.27], [0.18, 0.27], [-0.18, 0.27]]] });
  near(b.A, p.A, 1e-12); near(b.Iy, p.Iy, 1e-12); near(b.Iz, p.Iz, 1e-12);
  // Bredt: J = 4 Am²/∮ds/t with Am = 0.38 × 0.57
  near(b.J, 4 * (0.38 * 0.57) ** 2 / (2 * 0.38 / 0.03 + 2 * 0.57 / 0.02), 1e-12);
  const q = rectHollowSection({ h: 0.3, b: 0.2, t: 0.01 });
  near(q.A, 0.3 * 0.2 - 0.28 * 0.18, 1e-12);
});

test('rolled I-sections reproduce published table values (IPE 300, HEB 300)', () => {
  const ipe = iSectionByName('IPE 300');
  near(ipe.A, 53.81e-4, 2e-3); near(ipe.Iy, 8356e-8, 2e-3); near(ipe.Iz, 603.8e-8, 2e-3);
  near(ipe.Wely, 557.1e-6, 2e-3); near(ipe.Wply, 628.4e-6, 2e-3); near(ipe.Wplz, 125.2e-6, 3e-3);
  near(ipe.Asz, 25.68e-4, 3e-3, 'shear area Avz');
  near(ipe.Iw, 125.9e-9, 5e-3, 'warping constant');
  near(ipe.iy, 0.1246, 2e-3); near(ipe.iz, 0.0335, 2e-3);
  const heb = iSectionByName('HEB 300');
  near(heb.A, 149.1e-4, 2e-3); near(heb.Iy, 25170e-8, 2e-3); near(heb.Iz, 8563e-8, 2e-3);
  near(heb.Wply, 1869e-6, 2e-3); near(heb.Asz, 47.43e-4, 3e-3);
  // areas of the remaining library profiles [cm²]
  const areas = { 'IPE 200': 28.48, 'IPE 400': 84.46, 'IPE 500': 115.5, 'IPE 600': 156.0, 'HEA 200': 53.83, 'HEA 300': 112.5, 'HEB 200': 78.08, 'HEB 400': 197.8, 'HEB 600': 270.0 };
  for (const [name, a] of Object.entries(areas)) near(iSectionByName(name).A, a * 1e-4, 3e-3, name);
  assert.equal(Object.keys(EURO_I_PROFILES).length, 11);
  assert.equal(rolledISection({ h: 0.3, b: 0.15, tw: 0.007, tf: 0.011 }).r, 0.007);
});

test('welded I, plate girder, tee and stiffened plate vs polygon integration', () => {
  const w = weldedISection({ bft: 0.3, tft: 0.02, hw: 0.76, tw: 0.012, bfb: 0.4, tfb: 0.03 });
  const poly = polygonSection({
    outer: [[-0.2, 0], [0.2, 0], [0.2, 0.03], [0.006, 0.03], [0.006, 0.79], [0.15, 0.79], [0.15, 0.81], [-0.15, 0.81], [-0.15, 0.79], [-0.006, 0.79], [-0.006, 0.03], [-0.2, 0.03]],
  });
  near(w.A, poly.A, 1e-12); near(w.zc, poly.zc, 1e-12); near(w.Iy, poly.Iy, 1e-12); near(w.Iz, poly.Iz, 1e-12);
  // plastic neutral axis divides the area in two: bottom flange 0.012 m², web 0.00912, top 0.006 → half = 0.01356 → 0.03 + 0.00156/0.012
  near(w.zp, 0.03 + 0.00156 / 0.012, 1e-9);
  // doubly symmetric case agrees with iSection
  const s = weldedISection({ bft: 0.2, tft: 0.015, hw: 0.37, tw: 0.008 });
  const i = iSection({ h: 0.4, b: 0.2, tw: 0.008, tf: 0.015 });
  near(s.Wply, i.Wply, 1e-9); near(s.Iy, i.Iy, 1e-12); near(s.Iw, i.Iw, 2e-3, 'flange-only vs full Iz in the warping constant'); near(s.zsc, 0.2, 1e-12);
  assert.equal(plateGirderSection({ bft: 0.3, tft: 0.02, hw: 1.2, tw: 0.01, stiffenerSpacing: 1.8 }).panelAspect, 1.5);
  const t = teeSection({ h: 0.2, b: 0.15, tw: 0.01, tf: 0.015 });
  const tp = polygonSection({ outer: [[-0.005, 0], [0.005, 0], [0.005, 0.185], [0.075, 0.185], [0.075, 0.2], [-0.075, 0.2], [-0.075, 0.185], [-0.005, 0.185]] });
  near(t.A, tp.A, 1e-12); near(t.zc, tp.zc, 1e-12); near(t.Iy, tp.Iy, 1e-12);
  const sp = stiffenedPlateStrip({ s: 0.7, tp: 0.016, hw: 0.3, tw: 0.012, bf: 0.1, tf: 0.016 });
  // hand: areas 0.0112, 0.0036, 0.0016 at z = 0.008, 0.166, 0.324 → zc = (8.96e-5 + 5.976e-4 + 5.184e-4)/0.0164
  near(sp.zc, (0.0112 * 0.008 + 0.0036 * 0.166 + 0.0016 * 0.324) / 0.0164, 1e-12);
  near(sp.WelPlate * sp.zc, sp.Iy, 1e-12);
  assert.ok(sp.WelStiff < sp.WelPlate);
});

test('channel and angle', () => {
  const c = channelSection({ h: 0.3, b: 0.1, tw: 0.01, tf: 0.016 });
  const p = polygonSection({ outer: [[0, -0.15], [0.1, -0.15], [0.1, -0.134], [0.01, -0.134], [0.01, 0.134], [0.1, 0.134], [0.1, 0.15], [0, 0.15]] });
  near(c.A, p.A, 1e-12); near(c.Iy, p.Iy, 1e-12); near(c.Iz, p.Iz, 1e-12); near(c.yc, p.yc, 1e-12);
  // shear centre e = 3 b'² tf/(6 b' tf + h' tw), b' = 0.095, h' = 0.284
  near(c.eShear, 3 * 0.095 ** 2 * 0.016 / (6 * 0.095 * 0.016 + 0.284 * 0.01), 1e-12);
  const a = angleSection({ h: 0.1, b: 0.1, t: 0.01 });
  near(a.A, 1900e-6, 1e-12);
  near(a.yc, (1000 * 5 + 900 * 55) / 1900 * 1e-3, 1e-12);
  near(a.zc, a.yc, 1e-12);
  near(Math.abs(a.alpha), Math.PI / 4, 1e-9, 'equal angle: principal axes at 45°');
  near(a.Iu + a.Iv, a.Iy + a.Iz, 1e-12);
  assert.ok(a.Iv < a.Iy && a.Iy < a.Iu);
});

test('polygon section: triangle and orientation independence', () => {
  // right triangle b × h: A = bh/2, centroid (b/3, h/3), I about centroidal axis parallel to the base = b h³/36
  const tri = polygonSection({ outer: [[0, 0], [0.3, 0], [0, 0.6]] });
  near(tri.A, 0.09, 1e-12); near(tri.yc, 0.1, 1e-12); near(tri.zc, 0.2, 1e-12);
  near(tri.Iy, 0.3 * 0.6 ** 3 / 36, 1e-12); near(tri.Iz, 0.6 * 0.3 ** 3 / 36, 1e-12);
  near(tri.Iyz, -(0.3 ** 2) * 0.6 ** 2 / 72, 1e-12);
  const cw = polygonSection({ outer: [[0, 0.6], [0.3, 0], [0, 0]] });
  near(cw.Iy, tri.Iy, 1e-12);
  const flat = polygonSection({ outer: Float64Array.from([0, 0, 0.3, 0, 0, 0.6]) });
  near(flat.A, 0.09, 1e-12);
});

test('concrete-filled tube and tubular library', () => {
  const s = compositeTubeSection({ D: 1.0, t: 0.025, Es: 2.1e11, Ec: 3.0e10 });
  const st = tubeSection({ D: 1.0, t: 0.025 });
  const Ac = Math.PI * 0.95 ** 2 / 4; const Ic = Math.PI * 0.95 ** 4 / 64;
  near(s.EA, 2.1e11 * st.A + 3.0e10 * Ac, 1e-12);
  near(s.EI, 2.1e11 * st.Iy + 3.0e10 * Ic, 1e-12);
  near(s.massPerLength, 7850 * st.A + 2400 * Ac, 1e-12);
  near(compositeTubeSection({ D: 1.0, t: 0.025, Es: 2.1e11, Ec: 3.0e10, concreteFactor: 0.6 }).EI, 2.1e11 * st.Iy + 0.6 * 3.0e10 * Ic, 1e-12);
  // pile grouted inside a leg: inner tube steel adds, annulus is grout
  const g = compositeTubeSection({ D: 1.2, t: 0.03, Es: 2.1e11, Ec: 2.0e10, innerD: 1.0, innerT: 0.04 });
  near(g.As, tubeSection({ D: 1.2, t: 0.03 }).A + tubeSection({ D: 1.0, t: 0.04 }).A, 1e-12);
  near(g.Ac, Math.PI / 4 * (1.14 ** 2 - 1.0 ** 2), 1e-12);
  assert.ok(OFFSHORE_TUBULARS.length > 20);
  const t = tubularByName('914x25.4');
  near(t.D, 0.9144, 1e-12);
  assert.equal(resolveSection({ type: 'tube', D: 0.5, t: 0.02 }).type, 'tube');
  assert.throws(() => resolveSection({ type: 'hexagon' }), RangeError);
  assert.throws(() => tubularByName('1x1'), RangeError);
});
