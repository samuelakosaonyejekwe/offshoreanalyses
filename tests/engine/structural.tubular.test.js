import test from 'node:test';
import assert from 'node:assert/strict';
import {
  isoTubularStrengths, isoTubularCheck, apiLrfdTubularCheck, apiWsdTubularCheck, hoopBucklingCoefficient, tubularMemberCheck,
} from '../../src/engine/structural/tubular.js';
import { tubeSection } from '../../src/engine/structural/sections.js';

const near = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);
const E = 2.05e11;
const uc = (r, name) => r.checks.find((c) => c.name.startsWith(name)).uc;

test('ISO 19902: stocky member reaches yield in compression and the plastic moment in bending', () => {
  // D/t = 20, fy = 355 MPa: fxe = 0.6 E t/D = 6150 MPa → fy/fxe = 0.058 ≤ 0.170 → fyc = fy; fy D/(E t) = 0.0346 ≤ 0.0517 → fb = (Zp/Ze) fy
  const D = 0.5; const t = 0.025; const fy = 355e6;
  const s = isoTubularStrengths({ D, t, fy, E, L: 0.5 });
  const g = tubeSection({ D, t });
  near(s.fyc, fy, 1e-12);
  assert.ok(s.lambda < 0.05);
  near(s.fc, fy, 1e-3, 'column strength → yield');
  near(s.fb, g.Wply / g.Wely * fy, 1e-12);
  near(s.Mp, g.Wply * fy, 1e-12);
  near(s.fv, fy / Math.sqrt(3), 1e-12);
  // pure axial compression at the design resistance: UC = γR,c σc/fc = 1
  const N = -s.A * s.fc / 1.18;
  near(isoTubularCheck({ D, t, fy, E, L: 0.5, N }).governing.uc, 1, 1e-9);
  // pure tension: UC = 1.05 σt/fy
  near(uc(isoTubularCheck({ D, t, fy, E, L: 0.5, N: 0.5 * s.A * fy }), 'Axial tension'), 1.05 * 0.5, 1e-12);
});

test('ISO 19902: slender member follows 0.9 × Euler; column curve is continuous at λ = 1.34', () => {
  const D = 0.3; const t = 0.012; const fy = 355e6; const L = 25;
  const s = isoTubularStrengths({ D, t, fy, E, L });
  assert.ok(s.lambda > 1.34, `λ = ${s.lambda}`);
  near(s.Pc, 0.9 * Math.PI ** 2 * E * s.I / L ** 2, 1e-12, 'Pc = 0.9 π²EI/L²');
  near(s.fe, Math.PI ** 2 * E / (L / s.r) ** 2, 1e-12);
  // effective length factor enters as (K L)
  near(isoTubularStrengths({ D, t, fy, E, L, K: 0.5 }).lambda, s.lambda / 2, 1e-12);
  // continuity: 1 − 0.278·1.34² = 0.5008 vs 0.9/1.34² = 0.5012
  const Lc = (dl) => (1.34 + dl) * Math.PI * s.r / Math.sqrt(fy / E);
  const lo = isoTubularStrengths({ D, t, fy, E, L: Lc(-1e-6) }).fc;
  const hi = isoTubularStrengths({ D, t, fy, E, L: Lc(1e-6) }).fc;
  near(lo, 0.5008 * fy, 1e-3); near(hi, lo, 1.5e-3);
});

test('ISO 19902: local buckling and bending strength ranges (hand values)', () => {
  // D/t = 100, fy = 355: fxe = 0.6·205000/100 = 1230 MPa; fy/fxe = 0.28862 → fyc = (1.047 − 0.274·0.28862) fy = 0.96792 fy
  const s = isoTubularStrengths({ D: 2.0, t: 0.02, fy: 355e6, E });
  near(s.fxe, 1230e6, 1e-12);
  near(s.fyc, (1.047 - 0.274 * 355 / 1230) * 355e6, 1e-12);
  near(s.fyc / 355e6, 0.96792, 1e-4);
  // fy D/(E t) = 0.17317 > 0.1034 → fb = (0.94 − 0.76·0.17317)(Zp/Ze) fy
  const x = 355 * 100 / 205000;
  near(s.fb, (0.94 - 0.76 * x) * (s.Zp / s.Ze) * 355e6, 1e-12);
  // intermediate range, D/t = 40: x = 0.06927 → 1.13 − 2.58 x = 0.95129
  const m = isoTubularStrengths({ D: 0.8, t: 0.02, fy: 355e6, E });
  near(m.fb / (m.Zp / m.Ze * 355e6), 1.13 - 2.58 * 355 * 40 / 205000, 1e-12);
  // the two bending expressions meet at 0.0517 and 0.1034 (to the rounding of the published coefficients)
  near(1.13 - 2.58 * 0.0517, 1.0, 4e-3);
  near(1.13 - 2.58 * 0.1034, 0.94 - 0.76 * 0.1034, 3e-3);
  assert.ok(isoTubularStrengths({ D: 2.6, t: 0.02, fy: 355e6, E }).warnings.length > 0, 'D/t = 130 is flagged');
});

test('hoop buckling coefficient branches and hoop strength', () => {
  const D = 1.0; const t = 0.02; // D/t = 50
  // long cylinder: μ ≥ 1.6 D/t = 80 → Ch = 0.44 t/D = 0.0088
  near(hoopBucklingCoefficient({ D, t, L: 10 }).Ch, 0.0088, 1e-12);
  near(hoopBucklingCoefficient({ D, t, L: 10 }).mu, 100, 1e-12);
  // 0.825 D/t = 41.25 ≤ μ = 60 < 80 → 0.0088 + 0.21·50³/60⁴ = 0.0108255
  near(hoopBucklingCoefficient({ D, t, L: 6 }).Ch, 0.0088 + 0.21 * 125000 / 60 ** 4, 1e-12);
  // 1.5 ≤ μ = 20 < 41.25 → ISO 0.737/(20 − 0.579); API 0.736/(20 − 0.636)
  near(hoopBucklingCoefficient({ D, t, L: 2 }).Ch, 0.737 / 19.421, 1e-12);
  near(hoopBucklingCoefficient({ D, t, L: 2, code: 'API' }).Ch, 0.736 / 19.364, 1e-12);
  near(hoopBucklingCoefficient({ D, t, L: 0.25, code: 'API' }).Ch, 0.755 / (2.5 - 0.559), 1e-12);
  near(hoopBucklingCoefficient({ D, t, L: 0.1 }).Ch, 0.8, 1e-12);
  // fhe = 2 Ch E t/D; ring spacing 2 m: fhe = 2·0.037949·205000·0.02 = 311.2 MPa → 0.55 fy < fhe ≤ 2.44 fy → fh = 0.7 fy (fhe/fy)^0.4
  const s = isoTubularStrengths({ D, t, fy: 355e6, E, L: 10, Lr: 2 });
  near(s.fhe, 2 * (0.737 / 19.421) * E * 0.02, 1e-12);
  near(s.fh, 0.7 * 355e6 * (s.fhe / 355e6) ** 0.4, 1e-12);
  // hoop check: σh = p D/(2t), UC = 1.25 σh/fh
  const p = 2e6;
  near(uc(isoTubularCheck({ D, t, fy: 355e6, E, L: 10, Lr: 2, p }), 'Hoop'), 1.25 * (p * D / (2 * t)) / s.fh, 1e-12);
});

test('ISO 19902 combined checks without pressure (hand calculation)', () => {
  const D = 0.6; const t = 0.02; const fy = 355e6; const L = 12;
  const s = isoTubularStrengths({ D, t, fy, E, L });
  const N = -1.5e6; const My = 2.0e5; const Mz = 1.0e5;
  const sc = 1.5e6 / s.A; const sby = My / s.Ze; const sbz = Mz / s.Ze;
  const r = isoTubularCheck({ D, t, fy, E, L, N, My, Mz, Cmy: 0.85, Cmz: 0.85 });
  const stab = 1.18 * sc / s.fc + 1.05 / s.fb * Math.hypot(0.85 * sby / (1 - sc / s.fe), 0.85 * sbz / (1 - sc / s.fe));
  const str = 1.18 * sc / s.fyc + 1.05 * Math.hypot(sby, sbz) / s.fb;
  near(uc(r, 'Compression + bending (beam-column)'), stab, 1e-12);
  near(uc(r, 'Compression + bending (local strength)'), str, 1e-12);
  near(r.governing.uc, Math.max(stab, str), 1e-12);
  const rt = isoTubularCheck({ D, t, fy, E, L, N: 1.5e6, My, Mz });
  near(uc(rt, 'Tension + bending'), 1.05 * sc / fy + 1.05 * Math.hypot(sby, sbz) / s.fb, 1e-12);
  // shear and torsion: τ = 2V/A, τt = T D/(2 Ip)
  const rv = isoTubularCheck({ D, t, fy, E, L, V: 5e5, T: 1e5 });
  near(uc(rv, 'Beam shear'), 1.05 * (2 * 5e5 / s.A) / (fy / Math.sqrt(3)), 1e-12);
  near(uc(rv, 'Torsional shear'), 1.05 * (1e5 * D / (2 * s.Ip)) / (fy / Math.sqrt(3)), 1e-12);
  // custom partial factors
  near(uc(isoTubularCheck({ D, t, fy, E, L, N: 1e6, gamma: { t: 1.0 } }), 'Axial tension'), 1e6 / s.A / fy, 1e-12);
});

test('ISO 19902 with hydrostatic pressure: reduced strengths and limiting cases', () => {
  const D = 1.0; const t = 0.025; const fy = 355e6; const L = 15; const p = 1.5e6;
  const s = isoTubularStrengths({ D, t, fy, E, L });
  const sh = p * D / (2 * t); const sq = 0.5 * sh;
  const B = 1.25 * sh / s.fh; const eta = 5 - 4 * s.fh / fy;
  const red = Math.sqrt(1 + 0.09 * B * B - B ** (2 * eta)) - 0.3 * B;
  // tension + bending + pressure
  const N = 4e6; const My = 6e5;
  const r = isoTubularCheck({ D, t, fy, E, L, N, My, p });
  const st = N / s.A - sq;
  near(uc(r, 'Tension + bending + hydrostatic'), 1.05 * st / (fy * red) + 1.05 * (My / s.Ze) / (s.fb * red), 1e-12);
  // compression: fc,h formula and its continuity with 0.9 fyc/λ² at λ = 1.34/√(1 − 2σq/fyc)
  const a = 0.1;
  const lam = 1.34 / Math.sqrt(1 - 2 * a);
  const c0 = 1 - 0.278 * lam * lam;
  const fch = 0.5 * (c0 - 2 * a + Math.sqrt(c0 * c0 + 1.12 * lam * lam * a));
  near(fch, 0.9 / (lam * lam), 2e-3, 'fc,h is continuous at the transition slenderness');
  // zero pressure limit of the formula: fc,h → (1 − 0.278 λ²) fyc
  near(0.5 * (c0 + Math.sqrt(c0 * c0)), c0, 1e-12);
  const rc = isoTubularCheck({ D, t, fy, E, L, N: -3e6, My, p });
  const sc = 3e6 / s.A;
  near(uc(rc, 'Compression + bending + hydrostatic pressure (local'), 1.18 * (sc + sq) / s.fyc + 1.05 * (My / s.Ze) / (s.fb * red), 1e-12);
  const l2 = s.lambda * s.lambda; const aa = sq / s.fyc; const c1 = 1 - 0.278 * l2;
  const fchm = 0.5 * s.fyc * (c1 - 2 * aa + Math.sqrt(c1 * c1 + 1.12 * l2 * aa));
  near(uc(rc, 'Compression + bending + hydrostatic pressure (beam'), 1.18 * sc / fchm + 1.05 * (0.85 * My / s.Ze / (1 - sc / s.fe)) / (s.fb * red), 1e-12);
  // with capped-end forces already in N the same member state gives the same utilisation
  const rc2 = isoTubularCheck({ D, t, fy, E, L, N: -3e6 - sq * s.A, My, p, cappedEndIncluded: true });
  near(rc2.governing.uc, rc.governing.uc, 1e-12);
  assert.ok(rc.governing.uc > isoTubularCheck({ D, t, fy, E, L, N: -3e6, My }).governing.uc, 'pressure increases the utilisation');
});

test('API RP 2A-WSD allowable stresses (hand calculation) and interaction', () => {
  const D = 0.6; const t = 0.02; const fy = 345e6; const Ea = 2.0e11; const L = 10;
  const g = tubeSection({ D, t });
  const kl = L / g.iy;
  const Cc = Math.sqrt(2 * Math.PI ** 2 * Ea / fy);
  assert.ok(kl < Cc);
  const Fa = (1 - kl * kl / (2 * Cc * Cc)) * fy / (5 / 3 + 3 * kl / (8 * Cc) - kl ** 3 / (8 * Cc ** 3));
  const r = apiWsdTubularCheck({ D, t, fy, E: Ea, L, N: -1e6, My: 1.5e5 });
  near(r.strengths.Fa, Fa, 1e-12);
  // D/t = 30 is just above 10340/345 = 29.97 → second bending range
  near(r.strengths.Fb, (0.84 - 1.74 * 345 * 30 / 200000) * fy, 1e-12);
  const fa = 1e6 / g.A; const fb = 1.5e5 / g.Wely;
  const Fe = 12 * Math.PI ** 2 * Ea / (23 * kl * kl);
  near(uc(r, 'Compression + bending (stability)'), fa / Fa + 0.85 * fb / (1 - fa / Fe) / r.strengths.Fb, 1e-12);
  near(uc(r, 'Compression + bending (strength)'), fa / (0.6 * fy) + fb / r.strengths.Fb, 1e-12);
});

test('API RP 2A-WSD bending ranges, Euler range, tension and one-third increase', () => {
  const fy = 345e6; const Ea = 2.0e11;
  // D/t = 25 ≤ 10340/345 = 29.97 → 0.75 Fy
  near(apiWsdTubularCheck({ D: 0.5, t: 0.02, fy, E: Ea, L: 5 }).strengths.Fb, 0.75 * fy, 1e-12);
  // D/t = 40: [0.84 − 1.74 Fy D/(E t)] Fy = (0.84 − 0.12006) Fy
  near(apiWsdTubularCheck({ D: 0.8, t: 0.02, fy, E: Ea, L: 5 }).strengths.Fb, (0.84 - 1.74 * 345 * 40 / 200000) * fy, 1e-12);
  // D/t = 80 > 20680/345 = 59.9: [0.72 − 0.58 Fy D/(E t)] Fy
  near(apiWsdTubularCheck({ D: 1.6, t: 0.02, fy, E: Ea, L: 5 }).strengths.Fb, (0.72 - 0.58 * 345 * 80 / 200000) * fy, 1e-12);
  // slender column: Fa = 12 π² E/(23 (Kl/r)²)
  const g = tubeSection({ D: 0.3, t: 0.012 });
  const L = 20; const kl = L / g.iy;
  assert.ok(kl > Math.sqrt(2 * Math.PI ** 2 * Ea / fy));
  const r = apiWsdTubularCheck({ D: 0.3, t: 0.012, fy, E: Ea, L, N: -1e5 });
  near(r.strengths.Fa, 12 * Math.PI ** 2 * Ea / (23 * kl * kl), 1e-12);
  near(uc(r, 'Axial compression'), (1e5 / g.A) / r.strengths.Fa, 1e-12);
  // tension 0.6 Fy; extreme condition allowables × 4/3
  const rt = apiWsdTubularCheck({ D: 0.3, t: 0.012, fy, E: Ea, L, N: 1e6, allowableFactor: 4 / 3 });
  near(uc(rt, 'Axial tension'), (1e6 / g.A) / (0.6 * fy * 4 / 3), 1e-12);
  // low axial load uses Eq. 3.3.1-3
  const low = apiWsdTubularCheck({ D: 0.6, t: 0.02, fy, E: Ea, L: 6, N: -1e5, My: 2e5 });
  assert.ok(low.checks.some((c) => c.ref.includes('3.3.1-3')));
  // hydrostatic: Fhc branches, tension interaction A² + B² + 2ν|A|B
  const D = 1.0; const t = 0.025; const p = 1e6;
  const h = apiWsdTubularCheck({ D, t, fy, E: Ea, L: 15, N: 3e6, p });
  const gg = tubeSection({ D, t });
  const fh = p * D / (2 * t);
  const A = (3e6 / gg.A - 0.5 * fh) / fy * 1.67; const B = fh / h.strengths.Fhc * 2.0;
  near(uc(h, 'Tension + hydrostatic'), A * A + B * B + 0.6 * Math.abs(A) * B, 1e-12);
  near(h.strengths.Fhe, 2 * hoopBucklingCoefficient({ D, t, L: 15, code: 'API' }).Ch * Ea * t / D, 1e-12);
});

test('API RP 2A-LRFD column curve and cosine interaction; dispatcher', () => {
  const D = 0.6; const t = 0.02; const fy = 345e6; const Ea = 2.0e11; const L = 10;
  const g = tubeSection({ D, t });
  const lam = L / (Math.PI * g.iy) * Math.sqrt(fy / Ea);
  const Fcn = (1 - 0.25 * lam * lam) * fy;
  const N = -2e6; const My = 3e5;
  const r = apiLrfdTubularCheck({ D, t, fy, E: Ea, L, N, My });
  near(r.strengths.Fcn, Fcn, 1e-12);
  // fy D/(E t) = 0.05175 is just above 0.0517 → second bending range
  const fc = 2e6 / g.A; const fb = My / g.Wely; const Fbn = (1.13 - 2.58 * 345 * 30 / 200000) * g.Wply / g.Wely * fy;
  near(uc(r, 'Compression + bending (local'), 1 - Math.cos(Math.PI / 2 * fc / (0.85 * fy)) + fb / (0.95 * Fbn), 1e-12);
  const Fe = Math.PI ** 2 * Ea / (L / g.iy) ** 2;
  near(uc(r, 'Compression + bending (beam'), fc / (0.85 * Fcn) + 0.85 * fb / (1 - fc / (0.85 * Fe)) / (0.95 * Fbn), 1e-12);
  // Euler range λ ≥ √2: Fcn = Fy/λ²
  const s = apiLrfdTubularCheck({ D: 0.3, t: 0.012, fy, E: Ea, L: 25, N: -1e5 });
  near(s.strengths.Fcn, fy / s.strengths.lambda ** 2, 1e-12);
  const viaDispatcher = tubularMemberCheck({ code: 'API-LRFD', section: g, fy, L, forces: { N, My }, options: { E: Ea } });
  near(viaDispatcher.governing.uc, r.governing.uc, 1e-12);
  assert.equal(tubularMemberCheck({ section: g, fy, L, forces: { N } }).governing.ref.startsWith('ISO 19902'), true);
  assert.throws(() => tubularMemberCheck({ code: 'XYZ', section: g, fy, L, forces: {} }), RangeError);
});
