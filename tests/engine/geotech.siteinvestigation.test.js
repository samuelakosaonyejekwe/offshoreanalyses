import { test } from 'node:test';
import assert from 'node:assert/strict';
import { parseCptTable, fitStrengthLine, cptDesignValues, boreholeDesignValues, PRESSURE_UNITS, LENGTH_UNITS } from '../../src/engine/geotech/siteinvestigation.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} expected ${b}, got ${a}`);

test('CPT table: header names and units in brackets, comma separated', () => {
  const t = parseCptTable('Depth [m],qc [MPa],fs [kPa],u2 [kPa]\n0.5,0.30,6.0,10\n1.0,0.33,6.6,20\n1.5,0.36,7.2,30\n');
  assert.deepEqual(t.z, [0.5, 1, 1.5]); assert.deepEqual(t.qc, [300000, 330000, 360000]); assert.deepEqual(t.fs, [6000, 6600, 7200]); assert.deepEqual(t.u2, [10000, 20000, 30000]);
  assert.equal(t.qcIsQt, false); assert.equal(t.skipped, 0); assert.deepEqual(t.units, { depth: 'm', qc: 'MPa', fs: 'kPa', u2: 'kPa' }); assert.deepEqual(t.warnings, []);
});
test('CPT table: semicolons with decimal commas, a units row, qt instead of qc, rows out of order', () => {
  const t = parseCptTable('z;qt;fs\nft;tsf;psi\n3,28;4;1\n1,64;3;0,5\n', { qc: 'kPa' });
  // 1.64 ft = 0.49987 m; 3 tsf = 287 281.5 Pa; 0.5 psi = 3 447.4 Pa
  close(t.z[0], 1.64 * 0.3048, 1e-12); close(t.z[1], 3.28 * 0.3048, 1e-12); close(t.qc[0], 3 * 95760.5, 1e-6); close(t.fs[0], 0.5 * 6894.757, 1e-9);
  assert.equal(t.qcIsQt, true); assert.equal(t.u2, null); assert.equal(t.units.qc, 'tsf');
});
test('CPT table: no header takes depth, cone resistance, sleeve friction in the units given', () => {
  const t = parseCptTable('0.5 300 6\n1.0 330 6.6\n\n1.0 999 9\nabc def ghi\n1.5 360 7.2', { qc: 'kPa', fs: 'kPa' });
  assert.deepEqual(t.z, [0.5, 1, 1.5]); assert.deepEqual(t.qc, [300000, 330000, 360000]); assert.match(t.format, /without header/);
});
test('CPT table: GEF file with void values', () => {
  const gef = '#GEFID= 1, 1, 0\n#COLUMN= 3\n#COLUMNINFO= 1, m, sondeerlengte, 1\n#COLUMNINFO= 2, MPa, conusweerstand, 2\n#COLUMNINFO= 3, MPa, plaatselijke wrijving, 3\n#COLUMNSEPARATOR= ;\n#RECORDSEPARATOR= !\n#COLUMNVOID= 2, -9999\n#EOH=\n0.5;0.3;0.006;!\n1.0;-9999;0.006;!\n1.5;0.35;0.007;!\n';
  const t = parseCptTable(gef);
  assert.equal(t.format, 'GEF'); assert.deepEqual(t.z, [0.5, 1.5]); assert.deepEqual(t.qc, [300000, 350000]); assert.deepEqual(t.fs, [6000, 7000]); assert.equal(t.skipped, 1); assert.match(t.warnings[0], /1 row was skipped/);
});
test('CPT table: SCPT group of an AGS 4 file', () => {
  const ags = '"GROUP","SCPT"\n"HEADING","LOCA_ID","SCPG_TESN","SCPT_DPTH","SCPT_RES","SCPT_FRES","SCPT_PWP2"\n"UNIT","","","m","MPa","MPa","MPa"\n"TYPE","ID","X","2DP","3DP","4DP","4DP"\n"DATA","CPT1","1","0.50","0.300","0.0060","0.0100"\n"DATA","CPT1","1","1.00","0.330","0.0066","0.0200"\n';
  const t = parseCptTable(ags);
  assert.equal(t.format, 'AGS'); assert.deepEqual(t.z, [0.5, 1]); assert.deepEqual(t.qc, [300000, 330000]); assert.deepEqual(t.fs, [6000, 6600]); assert.deepEqual(t.u2, [10000, 20000]);
});
test('CPT table: wrong units and unreadable input are reported', () => {
  // cone resistance typed in kPa but read as MPa: median 300 000 MPa is outside what soils give
  assert.match(parseCptTable('depth,qc,fs\n1,300,6\n2,330,7\n3,360,8').warnings.join(' '), /check the units/);
  assert.throws(() => parseCptTable('depth,fs\n1,6\n2,7'), /cone resistance/); assert.throws(() => parseCptTable(''), /no data rows/); assert.throws(() => parseCptTable('depth,qc,fs\n1,0.3,6'), /fewer than two/);
  assert.throws(() => parseCptTable('depth,qc,fs\n1,0.3,6\n2,0.4,7', { qc: 'furlong' }), /unknown qc unit/);
  assert.equal(PRESSURE_UNITS.MPa, 1e6); assert.equal(LENGTH_UNITS.ft, 0.3048);
});

test('strength line: least squares, hand calculation', () => {
  // su = 10, 12, 14 kPa at 1, 2, 3 m lie on su = 8 + 2 z exactly
  const a = fitStrengthLine([1, 2, 3], [10e3, 12e3, 14e3]);
  close(a.su0, 8000, 1e-9); close(a.suGrad, 2000, 1e-9); close(a.r2, 1, 1e-12); assert.equal(a.throughOrigin, false);
  // (1, 9), (2, 13), (3, 14) kPa: mean z 2, mean su 12; slope = [(-1)(-3) + 0 + (1)(2)]/2 = 2.5; intercept 12 − 5 = 7
  // residuals −0.5, 1, −0.5 → SSE 1.5; total 14 → R² = 1 − 1.5/14 = 0.89286
  const b = fitStrengthLine([1, 2, 3], [9e3, 13e3, 14e3]);
  close(b.su0, 7000, 1e-9); close(b.suGrad, 2500, 1e-9); close(b.r2, 0.892857, 1e-6);
  // (1, 1), (2, 5), (3, 9) kPa: free fit gives −3 at the mudline, so the line through zero is used: k = Σzs/Σz² = (1 + 10 + 27)/14 = 2.7143
  const c = fitStrengthLine([1, 2, 3], [1e3, 5e3, 9e3]);
  assert.equal(c.su0, 0); close(c.suGrad, 38000 / 14, 1e-9); assert.equal(c.throughOrigin, true);
  // strength falling with depth: no gradient, the mean is used
  const d = fitStrengthLine([1, 2, 3], [9e3, 8e3, 7e3]); assert.equal(d.suGrad, 0); close(d.su0, 8000, 1e-9);
  assert.throws(() => fitStrengthLine([1], [1]), /two points/);
});

test('CPT reduction: soft clay sounding, hand calculation of the first two readings', () => {
  // qc = 0.30, 0.33, 0.36 MPa and fs = 6.0, 6.6, 7.2 kPa at 1, 2, 3 m; no pore pressure, Nkt = 15.
  // Reading 1: Rf = 2 %; γ/γw = 0.27 log 2 + 0.36 log 3 + 1.236 = 0.08128 + 0.17176 + 1.236 = 1.48904 → γ = 9810 × 1.48904 = 14 607.5 N/m³
  //   σv0 = 14 607.5 Pa; su = (300 000 − 14 607.5)/15 = 19 026.2 Pa.
  // Reading 2: Rf = 2 %; γ/γw = 0.08128 + 0.36 log 3.3 + 1.236 = 1.50394 → γ = 14 753.7; σv0 = 14 607.5 + 14 753.7 = 29 361.2 Pa
  //   su = (330 000 − 29 361.2)/15 = 20 042.6 Pa.
  const d = cptDesignValues({ z: [1, 2, 3], qc: [3e5, 3.3e5, 3.6e5], fs: [6e3, 6.6e3, 7.2e3], Nkt: 15 });
  close(d.points.gamma[0], 14607.5, 0.5); close(d.points.sigmaV0[1], 29361.2, 0.5); close(d.points.su[1], 20042.6, 0.5);
  // this near the mudline the first reading plots as a sand mixture (Ic 2.48), the next two as silt mixtures: too few for a line
  assert.deepEqual(d.points.group, ['sand', 'silt', 'silt']); assert.equal(d.su0, null); assert.equal(d.phiDeg, null);
  close(d.share.silt, 2 / 3, 1e-12); close(d.share.sand, 1 / 3, 1e-12); assert.equal(d.depth, 3);
  // with neither a clay line nor a sand angle, the unit weight is the mean of all readings: (14 607.5 + 14 753.7 + 14 887.1)/3 − 10 052 = 4 697.4 N/m³
  close(d.gammaSub, 4697.4, 0.5);
});
test('CPT reduction: clay over sand gives a strength line, a cautious friction angle and the unit weight of the main soil', () => {
  const z = Array.from({ length: 40 }, (_, i) => 0.25 * (i + 1));
  const clay = (v) => v <= 6, qc = z.map((v) => (clay(v) ? 0.2e6 + 0.03e6 * v : 8e6 + 0.8e6 * v)), fs = z.map((v) => (clay(v) ? 4e3 + 500 * v : 40e3 + 3e3 * v)), u2 = z.map((v) => (clay(v) ? 15e3 * v : 10e3 * v));
  const d = cptDesignValues({ z, qc, fs, u2, a: 0.8, Nkt: 15 });
  // qt = qc + 0.2 u2 rises by 33 kPa/m in the clay, total stress by about 14.5 kPa/m: su gradient ≈ (33 − 14.5)/15 ≈ 1.2 kPa/m
  assert.ok(d.suGrad > 1100 && d.suGrad < 1350, `suGrad ${d.suGrad}`); assert.ok(d.su0 > 12e3 && d.su0 < 15e3, `su0 ${d.su0}`); assert.ok(d.suFit.r2 > 0.99);
  // su is inversely proportional to Nkt: 15/20 and 15/10 of the fitted values
  close(d.suRange.su0[0], d.su0 * 0.75, 1e-9); close(d.suRange.suGrad[1], d.suGrad * 1.5, 1e-9);
  const again = cptDesignValues({ z, qc, fs, u2, a: 0.8, Nkt: 10 }); close(again.su0, d.su0 * 1.5, 1e-6); close(again.suGrad, d.suGrad * 1.5, 1e-6);
  // dense sand below 6 m: the lower quartile lies between the weakest reading and the median
  assert.ok(d.phiRange[0] <= d.phiDeg && d.phiDeg <= d.phiRange[1] && d.phiRange[1] <= d.phiRange[2]); assert.ok(d.phiRange[2] <= 45 && d.phiRange[1] > 40);
  assert.ok(d.share.sand > 0.5); assert.equal(d.gammaSub, d.gammaSubSand); assert.ok(d.gammaSubSand > d.gammaSubClay && d.gammaSubSand > 6e3 && d.gammaSubSand < 11e3);
  // a table that already holds qt is not corrected again
  const qt = cptDesignValues({ z, qc, fs, u2, qcIsQt: true, Nkt: 15 }); assert.equal(qt.points.qt[3], qc[3]); assert.ok(d.points.qt[3] > qc[3]);
});

test('borehole layers: strength line, lowest sand friction angle, unit weight of the main soil', () => {
  // one clay layer 5 → 15 kPa over 0–5 m: su = 5 + 2 z exactly; clay is 5 of 8 m
  const d = boreholeDesignValues([{ top: 0, bottom: 5, soil: 'clay', suTop: 5e3, suBottom: 15e3, gammaSub: 6e3 }, { top: 5, bottom: 8, soil: 'sand', phiDeg: 34, gammaSub: 9.5e3 }]);
  close(d.su0, 5000, 1e-6); close(d.suGrad, 2000, 1e-6); assert.equal(d.phiDeg, 34); assert.equal(d.gammaSub, 6000); assert.deepEqual(d.share, { clay: 0.625, sand: 0.375 }); assert.equal(d.depth, 8);
  const s = boreholeDesignValues([{ top: 0, bottom: 4, soil: 'sand', phiDeg: 36, gammaSub: 10e3 }, { top: 4, bottom: 10, soil: 'sand', phiDeg: 32, gammaSub: 9e3 }]);
  assert.equal(s.su0, null); assert.equal(s.phiDeg, 32); close(s.gammaSub, (4 * 10e3 + 6 * 9e3) / 10, 1e-9);
  assert.throws(() => boreholeDesignValues([]), /at least one layer/); assert.throws(() => boreholeDesignValues([{ top: 2, bottom: 1, soil: 'clay', suTop: 1, suBottom: 1 }]), /below top/);
  assert.throws(() => boreholeDesignValues([{ top: 0, bottom: 1, soil: 'peat' }]), /'clay' or 'sand'/); assert.throws(() => boreholeDesignValues([{ top: 0, bottom: 1, soil: 'sand', phiDeg: 80 }]), /friction angle/);
});
