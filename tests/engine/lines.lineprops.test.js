import test from 'node:test';
import assert from 'node:assert/strict';
import * as lp from '../../src/engine/lines/lineprops.js';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);

test('chain break load BL = c d^2 (44 - 0.08 d): R4 100 mm = 9864 kN, R3 76 mm = 4884 kN (IACS W22 catalogue values)', () => {
  // 0.0274 * 100^2 * (44 - 8) = 9864 kN
  assert.ok(rel(lp.chainMBL({ d: 0.1, grade: 'R4' }), 9864e3) < 1e-12);
  // 0.0223 * 76^2 * (44 - 6.08) = 4884.3 kN
  assert.ok(rel(lp.chainMBL({ d: 0.076, grade: 'R3' }), 0.0223 * 5776 * 37.92 * 1e3) < 1e-12);
  assert.ok(rel(lp.chainMBL({ d: 0.076, grade: 'R3' }), 4884e3) < 1e-3);
  // grade ordering and ratios follow the constants
  const g = ['R3', 'R3S', 'R4', 'R4S', 'R5'].map((grade) => lp.chainMBL({ d: 0.12, grade }));
  for (let i = 1; i < g.length; i++) assert.ok(g[i] > g[i - 1]);
  assert.ok(rel(g[4] / g[0], 0.0320 / 0.0223) < 1e-12);
  assert.throws(() => lp.chainMBL({ d: 0.1, grade: 'R6' }), /grade/);
  assert.throws(() => lp.chainMBL({ d: -1 }), /d/);
});

test('chain mass, stiffness and submerged weight (hand calculation, 100 mm studless / stud-link)', () => {
  const a = lp.chainProperties({ d: 0.1, grade: 'R4', type: 'studless' });
  assert.ok(rel(a.mass, 199) < 1e-12);       // 19.9e3 * 0.1^2 kg/m
  assert.ok(rel(a.EA, 0.854e9) < 1e-12);     // 0.854e11 * 0.1^2 N
  assert.ok(rel(a.w, 199 * 9.80665 * (1 - 1025 / 7850)) < 1e-12); // 1696.7 N/m
  // volume-equivalent diameter ~ 1.8 d
  assert.ok(Math.abs(a.dVol / 0.1 - 1.8) < 0.01);
  const b = lp.chainProperties({ d: 0.1, type: 'studlink' });
  assert.ok(rel(b.mass, 219) < 1e-12 && rel(b.EA, 1.01e9) < 1e-12);
  assert.equal(b.MBL, a.MBL);
  // overrides
  assert.equal(lp.chainProperties({ d: 0.1, massCoeff: 2e4 }).mass, 200);
  assert.equal(a.indicative, true);
});

test('wire rope: six-strand break load coincides with EN 12385-4 K d^2 Rr at grade 1770', () => {
  const iwrc = lp.wireRopeProperties({ d: 0.09 });
  assert.ok(rel(iwrc.MBL, 0.356 * 90 * 90 * 1770) < 0.01); // K = 0.356, d in mm, Rr in N/mm2 -> N
  const fc = lp.wireRopeProperties({ d: 0.09, construction: 'six-strand-fc' });
  assert.ok(rel(fc.MBL, 0.330 * 90 * 90 * 1770) < 0.01);
  assert.ok(iwrc.mass > fc.mass && iwrc.EA > fc.EA);
  const sp = lp.wireRopeProperties({ d: 0.09, construction: 'spiral-strand' });
  assert.ok(sp.MBL > iwrc.MBL && sp.EA > 2 * iwrc.EA);
  // metallic fill factor implied by the mass coefficient must be physical (0.5 - 0.8)
  for (const c of Object.values(lp.WIRE_COEFFS)) {
    const fill = c.massCoeff / (7850 * Math.PI / 4);
    assert.ok(fill > 0.5 && fill < 0.8);
  }
  // sheathed rope: buoyancy on the sheath envelope
  const sh = lp.wireRopeProperties({ d: 0.09, construction: 'spiral-strand', dOuter: 0.106 });
  assert.ok(rel(sh.w, (sp.mass - 1025 * Math.PI * 0.106 ** 2 / 4) * 9.80665) < 1e-12);
});

test('fibre rope properties and load-dependent stiffness EA = (alpha + beta Tmean/MBL) MBL', () => {
  const pe = lp.fibreRopeProperties({ d: 0.2 });
  assert.ok(rel(pe.mass, 0.7978e3 * 0.04) < 1e-12);
  assert.ok(rel(pe.w, pe.mass * 9.80665 * (1 - 1025 / 1380)) < 1e-12);
  assert.ok(lp.fibreRopeProperties({ d: 0.2, material: 'hmpe' }).w < 0); // HMPE floats
  // polyester dynamic stiffness at 20 % MBL mean load: Kr = 18.5 + 0.33 * 20 = 25.1
  const k = lp.fibreRopeStiffness({ MBL: 1e7, Tmean: 2e6, preset: 'polyester-dynamic' });
  assert.ok(rel(k.Kr, 25.1) < 1e-12 && rel(k.EA, 2.51e8) < 1e-12);
  assert.equal(lp.fibreRopeStiffness({ MBL: 1e7, Tmean: 2e6, alpha: 10, beta: 20 }).Kr, 14);
  assert.throws(() => lp.fibreRopeStiffness({ MBL: 1e7, Tmean: 2e6 }), /alpha/);
});

test('corrosion allowance, DNV-OS-E301 drag coefficients and marine growth (hand calculation)', () => {
  // 0.4 mm/yr over 25 years on 120 mm chain -> 110 mm
  const d = lp.corrodedDiameter({ d: 0.12, rate: 0.4e-3, years: 25 });
  assert.ok(rel(d, 0.11) < 1e-12);
  assert.ok(lp.chainMBL({ d }) < lp.chainMBL({ d: 0.12 }));
  assert.deepEqual({ ...lp.HYDRO_COEFFS['studless-chain'] }, { Cdn: 2.4, Cdt: 1.15 });
  assert.deepEqual({ ...lp.HYDRO_COEFFS['stud-chain'] }, { Cdn: 2.6, Cdt: 1.4 });
  assert.equal(lp.hydroCoefficients({ type: 'six-strand-rope' }).Cdn, 1.8);
  assert.equal(lp.hydroCoefficients({ type: 'spiral-rope-sheathed' }).Cdn, 1.2);
  assert.equal(lp.hydroCoefficients({ type: 'spiral-rope-sheathed', ropeCdt: 0.02 }).Cdt, 0.02);
  // marine growth: D = 0.1 m, t = 0.05 m, rho = 1325, chain (mu = 2):
  // M = pi/4 (0.2^2 - 0.1^2) * 1325 * 2 = 62.44 kg/m; W = M (1 - 1025/1325) g = 138.6 N/m; Cd factor 2
  const mg = lp.marineGrowth({ d: 0.1, thickness: 0.05, Cdn: 2.4 });
  assert.ok(rel(mg.mass, Math.PI / 4 * 0.03 * 1325 * 2) < 1e-12 && rel(mg.mass, 62.44) < 1e-3);
  assert.ok(rel(mg.w, 62.439 * (300 / 1325) * 9.80665) < 1e-4);
  assert.ok(rel(mg.CdnGrowth, 4.8) < 1e-12);
  assert.ok(rel(lp.marineGrowth({ d: 0.1, thickness: 0.05, component: 'rope' }).mass, mg.mass / 2) < 1e-12);
  const seg = lp.toSegment(lp.chainProperties({ d: 0.1 }), 500, 'chain');
  assert.equal(seg.L, 500); assert.ok(seg.w > 0 && seg.EA > 0);
});
