import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  STEFAN_BOLTZMANN, isoFireC, hydrocarbonFireC, externalFireC, fireCurve, steelSpecificHeat, steelConductivity, steelThermalStrain,
  steelReductionFactors, criticalTemperatureC, tubeSectionFactor, steelTemperatureUnprotected, steelTemperatureProtected,
  timeToTemperature, pfpThickness, steelAtTemperature, restrainedThermalStress, steadyMultilayer, lumpedCooling, conduction1d,
  semiInfiniteStep,
} from '../../src/engine/structural/thermal.js';
import { erfc } from '../../src/engine/math/special.js';
import { rk45 } from '../../src/engine/math/ode.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b} (rel ${rel})`);
const near = (a, b, abs, msg) => assert.ok(Math.abs(a - b) <= abs, `${msg || ''} got ${a}, expected ${b} (abs ${abs})`);

test('ISO 834 standard fire: 842 °C at 30 min, 945 °C at 60 min, 1049 °C at 120 min (EN 1991-1-2 Eq. 3.4)', () => {
  near(isoFireC(1800), 842, 0.5); near(isoFireC(3600), 945, 0.5); near(isoFireC(7200), 1049, 0.5);
  assert.equal(isoFireC(0), 20);
  near(isoFireC(1800), 20 + 345 * Math.log10(241), 1e-9);
});

test('hydrocarbon and external fire curves: closed-form values and plateaus 1100 °C / 680 °C (EN 1991-1-2 Eq. 3.5, 3.6)', () => {
  near(hydrocarbonFireC(0), 20, 1e-9); near(externalFireC(0), 20, 1e-9);
  // t = 10 min: 20 + 1080(1 − 0.325 e^−1.67 − 0.675 e^−25) = 20 + 1080·(1 − 0.06119) = 1033.9 °C
  near(hydrocarbonFireC(600), 1033.9, 0.2);
  near(hydrocarbonFireC(4 * 3600), 1100, 1e-6);
  // t = 10 min: 20 + 660(1 − 0.687 e^−3.2) = 20 + 660·(1 − 0.02800) = 661.5 °C
  near(externalFireC(600), 661.5, 0.2);
  near(externalFireC(3 * 3600), 680, 1e-6);
  assert.equal(fireCurve('hydrocarbon', 0).alphaC, 50);
  assert.equal(fireCurve('iso834', 60).alphaC, 25);
  assert.throws(() => fireCurve('jet', 0), RangeError);
});

test('EN 1993-1-2 Table 3.1: ky,θ and kE,θ at the tabulated temperatures, interpolation between', () => {
  const ky = { 20: 1, 400: 1, 500: 0.78, 600: 0.47, 700: 0.23, 800: 0.11, 900: 0.06, 1000: 0.04, 1100: 0.02, 1200: 0 };
  const kE = { 100: 1, 200: 0.9, 300: 0.8, 400: 0.7, 500: 0.6, 600: 0.31, 700: 0.13, 800: 0.09, 900: 0.0675, 1000: 0.045, 1100: 0.0225, 1200: 0 };
  for (const [T, v] of Object.entries(ky)) near(steelReductionFactors(+T).ky, v, 1e-12, `ky(${T})`);
  for (const [T, v] of Object.entries(kE)) near(steelReductionFactors(+T).kE, v, 1e-12, `kE(${T})`);
  near(steelReductionFactors(550).ky, 0.625, 1e-12); near(steelReductionFactors(650).kE, 0.22, 1e-12);
  near(steelReductionFactors(300).kp, 0.613, 1e-12);
  assert.deepEqual(steelReductionFactors(1300), { ky: 0, kp: 0, kE: 0 });
  const s = steelAtTemperature({ fy: 355e6, thetaC: 600 });
  close(s.fy, 0.47 * 355e6, 1e-12); close(s.E, 0.31 * 2.1e11, 1e-12);
});

test('EN 1993-1-2 Eq. (4.22): critical temperature 585 °C at μ0 = 0.5, 526 °C at 0.7; table inverse of ky', () => {
  near(criticalTemperatureC(0.5), 585, 0.5);
  // μ0 = 0.7: 0.7^3.833 = 0.2548, 1/(0.9674·0.2548) − 1 = 3.057, ln = 1.1175 → 39.19·1.1175 + 482 = 525.8
  near(criticalTemperatureC(0.7), 525.8, 0.3);
  near(criticalTemperatureC(0.47, 'table'), 600, 1e-9);
  near(criticalTemperatureC(0.625, 'table'), 550, 1e-9);
  near(steelReductionFactors(criticalTemperatureC(0.3, 'table')).ky, 0.3, 1e-12);
  assert.throws(() => criticalTemperatureC(0.001), RangeError);
});

test('EN 1993-1-2 §3.4: carbon steel specific heat (440 at 20 °C, peak 5000 at 735 °C, 650 above 900 °C), conductivity, elongation', () => {
  near(steelSpecificHeat(20), 425 + 0.773 * 20 - 1.69e-3 * 400 + 2.22e-6 * 8000, 1e-9);
  near(steelSpecificHeat(20), 439.8, 0.05);
  near(steelSpecificHeat(735), 5000, 1e-6); // 545 + 17820/4 = 666 + 13002/3
  near(666 + 13002 / 3, 5000, 1e-9);
  near(steelSpecificHeat(600), 666 + 13002 / 138, 1e-9);
  // continuity at 600 °C of the two branches: 425 + 463.8 − 608.4 + 479.5 = 759.9 vs 760.2
  near(425 + 0.773 * 600 - 1.69e-3 * 36e4 + 2.22e-6 * 216e6, 666 + 13002 / 138, 0.5);
  assert.equal(steelSpecificHeat(1000), 650);
  near(steelConductivity(20), 53.334, 1e-9); assert.equal(steelConductivity(900), 27.3);
  near(steelThermalStrain(20), 0, 2e-6); near(steelThermalStrain(800), 0.011, 1e-12);
  near(steelThermalStrain(600), 1.2e-5 * 600 + 0.4e-8 * 36e4 - 2.416e-4, 1e-15);
  near(steelThermalStrain(1000), 0.0138, 1e-12);
});

test('unprotected steel, EN 1993-1-2 §4.2.5.1: step equation agrees with an adaptive ODE solution; bounded by the gas temperature', () => {
  const AmV = 150;
  const r = steelTemperatureUnprotected({ sectionFactor: AmV, fire: 'iso834', duration: 1800, dt: 1 });
  // independent integration of dθ/dt = (Am/V)/(ca ρa) [αc(θg − θ) + εm σ((θg + 273)⁴ − (θ + 273)⁴)]
  const f = (t, y) => {
    const g = isoFireC(t);
    return [AmV / (steelSpecificHeat(y[0]) * 7850) * (25 * (g - y[0]) + 0.7 * STEFAN_BOLTZMANN * ((g + 273) ** 4 - (y[0] + 273) ** 4))];
  };
  const ode = rk45(f, 0, [20], 1800, { rtol: 1e-8, atol: 1e-8, hMax: 5 });
  near(r.finalC, ode.y[ode.y.length - 1][0], 1.5, 'step equation vs ODE at 30 min');
  assert.ok(r.finalC < isoFireC(1800) && r.finalC > 700, `θa(30 min) = ${r.finalC}`);
  for (let i = 0; i < r.t.length; i++) assert.ok(r.steelC[i] <= r.gasC[i] + 1e-9);
  // heavier section (smaller Am/V) heats more slowly; shadow factor scales the section factor
  const heavy = steelTemperatureUnprotected({ sectionFactor: 50, fire: 'iso834', duration: 1800 });
  assert.ok(heavy.finalC < r.finalC);
  const sh = steelTemperatureUnprotected({ sectionFactor: 200, ksh: 0.75, fire: 'iso834', duration: 1800, dt: 1 });
  near(sh.finalC, r.finalC, 1e-9);
  // hydrocarbon fire is more severe (αc = 50 and faster rise)
  assert.ok(steelTemperatureUnprotected({ sectionFactor: AmV, fire: 'hydrocarbon', duration: 600 }).finalC
    > steelTemperatureUnprotected({ sectionFactor: AmV, fire: 'iso834', duration: 600 }).finalC);
  // time to a critical temperature is consistent with the history
  const tc = timeToTemperature(r, 585);
  assert.ok(tc > 0 && tc < 1800);
  near(steelTemperatureUnprotected({ sectionFactor: AmV, fire: 'iso834', duration: tc, dt: 1 }).finalC, 585, 1.5);
  assert.equal(timeToTemperature(r, 2000), Infinity);
  assert.throws(() => steelTemperatureUnprotected({ sectionFactor: AmV, fire: 'iso834', duration: 600, dt: 10 }), RangeError);
  assert.throws(() => steelTemperatureUnprotected({ sectionFactor: AmV, fire: (t) => 500, duration: 600 }), RangeError);
  close(tubeSectionFactor({ D: 0.5, t: 0.02 }), 0.5 / (0.48 * 0.02), 1e-12); // ≈ 1/t for thin walls
});

test('insulated steel, EN 1993-1-2 §4.2.5.2: light insulation (φ = 0) in a constant-temperature furnace follows the exponential solution', () => {
  // φ = 0, constant ca would give θ = θg − (θg − θ0) exp(−λp Ap/V t/(dp ca ρa)); with ca(θ) compare against the ODE
  const o = { sectionFactor: 120, dp: 0.02, lambdaP: 0.15, fire: () => 800, duration: 3600, dt: 5 };
  const r = steelTemperatureProtected(o);
  const ode = rk45((t, y) => [0.15 * 120 * (800 - y[0]) / (0.02 * steelSpecificHeat(y[0]) * 7850)], 0, [20], 3600, { rtol: 1e-9, atol: 1e-9 });
  near(r.finalC, ode.y[ode.y.length - 1][0], 1.0);
  // heat capacity of the protection delays heating; thicker protection lowers the temperature
  const heavy = steelTemperatureProtected({ ...o, fire: 'iso834', rhoP: 800, cP: 1100, dt: 30 });
  const light = steelTemperatureProtected({ ...o, fire: 'iso834', dt: 30 });
  assert.ok(heavy.finalC < light.finalC);
  assert.ok(steelTemperatureProtected({ ...o, fire: 'iso834', dp: 0.04 }).finalC < light.finalC);
  // never negative increments while the gas temperature rises
  for (let i = 1; i < heavy.t.length; i++) assert.ok(heavy.steelC[i] >= heavy.steelC[i - 1] - 1e-12);
  assert.throws(() => steelTemperatureProtected({ ...o, dt: 60 }), RangeError);
});

test('passive fire protection sizing: thickness found by root finding reproduces the limiting temperature', () => {
  const o = { sectionFactor: 150, lambdaP: 0.12, rhoP: 600, cP: 1000, fire: 'hydrocarbon', duration: 7200, limitC: 400 };
  const r = pfpThickness(o);
  assert.ok(r.required && r.dp > 0.005 && r.dp < 0.2, `dp = ${r.dp}`);
  near(r.steelC, 400, 0.01);
  // a longer rating or a lower limit needs more protection
  assert.ok(pfpThickness({ ...o, duration: 10800 }).dp > r.dp);
  assert.ok(pfpThickness({ ...o, limitC: 300 }).dp > r.dp);
  assert.throws(() => pfpThickness({ ...o, limitC: 25, dpMax: 0.01 }), RangeError);
});

test('restrained thermal stress: σ = −E α ΔT = −252 MPa for ΔT = 100 K in steel; spring restraint; yield cap', () => {
  const r = restrainedThermalStress({ E: 2.1e11, dT: 100 });
  close(r.stress, -252e6, 1e-12);
  // restraint as stiff as the member: half the stress
  const A = 0.01; const L = 5;
  const half = restrainedThermalStress({ E: 2.1e11, dT: 100, A, L, ks: 2.1e11 * A / L });
  close(half.stress, -126e6, 1e-12); close(half.force, -1.26e6, 1e-12);
  const y = restrainedThermalStress({ E: 2.1e11, dT: 300, fy: 355e6 });
  assert.ok(y.yielded); close(y.stress, -355e6, 1e-12); close(y.elasticStress, -756e6, 1e-12);
  // EN 1993-1-2 strain: 20 → 500 °C free strain 6.758e-3
  close(restrainedThermalStress({ E: 0.6 * 2.1e11, thermalStrain: steelThermalStrain(500) - steelThermalStrain(20) }).freeStrain, 6.76e-3, 1e-3);
});

test('steady multilayer wall and cylinder: series resistances, hand calculation', () => {
  // wall: 1/10 + 0.2/1.0 + 0.05/0.04 + 1/25 = 0.1 + 0.2 + 1.25 + 0.04 = 1.59 m²K/W; ΔT = 40 K → q = 25.157 W/m²
  const w = steadyMultilayer({ layers: [{ thickness: 0.2, k: 1.0 }, { thickness: 0.05, k: 0.04 }], hInner: 10, hOuter: 25, TinnerC: 20, TouterC: -20 });
  close(w.resistance, 1.59, 1e-12); close(w.q, 40 / 1.59, 1e-12); close(w.U, 1 / 1.59, 1e-12);
  close(w.interfaceC[0], 20 - 2.5157, 1e-4); close(w.interfaceC[2], -20 + 25.157 * 0.04, 1e-4);
  // pipe 0.1 m bore radius, 10 mm steel (k = 45), 40 mm insulation (k = 0.2), 60 °C inside wall, sea water h = 500 at 4 °C
  // R' = ln(0.11/0.1)/(2π·45) + ln(0.15/0.11)/(2π·0.2) + 1/(2π·0.15·500) = 3.371e-4 + 0.24682 + 2.122e-3 = 0.24928 mK/W
  const c = steadyMultilayer({ geometry: 'cylinder', rInner: 0.1, layers: [{ thickness: 0.01, k: 45 }, { thickness: 0.04, k: 0.2 }], hOuter: 500, TinnerC: 60, TouterC: 4 });
  close(c.resistance, 0.24928, 1e-4); close(c.q, 56 / 0.24928, 1e-4);
  close(c.U, 1 / (0.24928 * 2 * Math.PI * 0.15), 1e-4);
});

test('lumped cooling: T = T∞ + (T0 − T∞) e^(−t/τ), τ = m c/(hA); Biot number', () => {
  // steel sphere D = 20 mm: m = 7850·4.189e-6 = 0.03288 kg, A = 1.2566e-3 m², c = 460, h = 30 → τ = 401.2 s
  const o = { T0C: 500, TinfC: 20, h: 30, area: Math.PI * 0.02 ** 2, mass: 7850 * Math.PI / 6 * 0.02 ** 3, cp: 460, k: 45, Lc: 0.02 / 6 };
  const r = lumpedCooling({ ...o, t: 401.2 });
  close(r.tau, 7850 * 460 * 0.02 / 6 / 30, 1e-12); close(r.tau, 401.2, 1e-3);
  near(r.TC, 20 + 480 / Math.E, 0.05);
  close(r.biot, 30 * 0.02 / 6 / 45, 1e-12); assert.ok(r.valid);
  close(r.timeTo(100), -r.tau * Math.log(80 / 480), 1e-12);
  assert.throws(() => r.timeTo(10), RangeError);
});

test('transient conduction, plane: semi-infinite solid with a surface temperature step, T = Ts + (T0 − Ts) erf(x/2√(αt)) (Carslaw & Jaeger §2.4)', () => {
  const k = 1.6; const rho = 2300; const cp = 900; const alpha = k / (rho * cp);
  const t = 3600; // penetration ≈ 4√(αt) = 0.21 m ≪ 0.6 m slab
  const o = { layers: [{ thickness: 0.6, k, rho, cp, n: 300 }], T0C: 20, left: { type: 'temperature', TC: 500 }, tEnd: t, dt: 10, theta: 0.5 };
  const r = conduction1d(o);
  for (const xq of [0.01, 0.02, 0.05, 0.1]) {
    const i = Math.round(xq / 0.002);
    const exact = 500 - 480 * (1 - erfc(xq / (2 * Math.sqrt(alpha * t))));
    near(r.finalC[i], exact, 0.6, `x = ${xq}`);
    near(semiInfiniteStep({ x: xq, t, alpha, T0C: 20, TsC: 500 }).TC, exact, 1e-9);
  }
  near(r.finalC[r.finalC.length - 1], 20, 1e-6);
  // surface flux k ΔT/√(παt)
  close(semiInfiniteStep({ x: 0, t, alpha, T0C: 20, TsC: 500, k }).flux, 1.6 * 480 / Math.sqrt(Math.PI * alpha * t), 1e-12);
  // convergence: backward Euler is first order in time — halving dt halves the error at x = 20 mm
  const exact = 500 - 480 * (1 - erfc(0.02 / (2 * Math.sqrt(alpha * t))));
  const e1 = Math.abs(conduction1d({ ...o, theta: 1, dt: 120 }).finalC[10] - exact);
  const e2 = Math.abs(conduction1d({ ...o, theta: 1, dt: 60 }).finalC[10] - exact);
  assert.ok(e1 / e2 > 1.6 && e1 / e2 < 2.4, `error ratio ${e1 / e2}`);
});

test('transient conduction, plane: energy conservation under a prescribed flux; insulated slab with convection reaches the film temperature', () => {
  const layer = { thickness: 0.05, k: 45, rho: 7850, cp: 460, n: 20 };
  const r = conduction1d({ layers: [layer], T0C: 20, left: { type: 'flux', q: 20e3 }, tEnd: 600, dt: 5 });
  // stored energy rise = q t exactly (backward Euler conserves energy)
  close(r.energy[r.energy.length - 1] - r.energy[0], 20e3 * 600, 1e-10);
  close(r.qLeft[r.qLeft.length - 1], 20e3, 1e-12);
  // mean temperature rise = q t/(ρ c L) = 66.46 K
  let mean = 0;
  for (let i = 0; i < 20; i++) mean += 0.5 * (r.finalC[i] + r.finalC[i + 1]) / 20;
  near(mean - 20, 20e3 * 600 / (7850 * 460 * 0.05), 1e-6);
  // thin plate (Bi = hL/k = 0.028) with convection on one face ≈ lumped solution
  const c = conduction1d({ layers: [layer], T0C: 400, right: { type: 'convection', h: 25, TinfC: 20 }, tEnd: 7200, dt: 10, theta: 0.5 });
  const lump = lumpedCooling({ T0C: 400, TinfC: 20, h: 25, area: 1, mass: 7850 * 0.05, cp: 460, t: 7200 });
  near(c.finalC[10], lump.TC, 1.5);
});

test('transient conduction, cylinder: multilayer pipe converges to the steady resistance solution; solid cylinder conserves energy', () => {
  const layers = [{ thickness: 0.01, k: 45, rho: 7850, cp: 460, n: 4 }, { thickness: 0.04, k: 0.2, rho: 900, cp: 2000, n: 16 }];
  const r = conduction1d({
    geometry: 'cylinder', rInner: 0.1, layers, T0C: 4, left: { type: 'temperature', TC: 60 }, right: { type: 'convection', h: 500, TinfC: 4 },
    tEnd: 40 * 3600, dt: 120,
  });
  const s = steadyMultilayer({ geometry: 'cylinder', rInner: 0.1, layers, hOuter: 500, TinnerC: 60, TouterC: 4 });
  near(r.finalC[r.interfaces[1]], s.interfaceC[1], 1e-3, 'steel / insulation interface');
  near(r.finalC[r.finalC.length - 1], s.interfaceC[2], 1e-3, 'outer surface');
  // outer surface flux × 2π r = heat loss per unit length
  close(-r.qRight[r.qRight.length - 1] * 2 * Math.PI * 0.15, s.q, 1e-4);
  // solid cylinder (r = 0 axis) heated by a surface flux: energy per unit length = q·2πR·t
  const sc = conduction1d({ geometry: 'cylinder', layers: [{ thickness: 0.05, k: 45, rho: 7850, cp: 460, n: 25 }], T0C: 20, right: { type: 'flux', q: 10e3 }, tEnd: 300, dt: 2 });
  close(sc.energy[sc.energy.length - 1] - sc.energy[0], 10e3 * 2 * Math.PI * 0.05 * 300, 1e-10);
  assert.ok(sc.finalC[0] < sc.finalC[25]);
  assert.throws(() => conduction1d({ geometry: 'cylinder', layers, T0C: 4, left: { type: 'temperature', TC: 60 }, tEnd: 10, dt: 1 }), RangeError);
});

test('transient conduction with radiation: plate in a furnace approaches the furnace temperature; thin plate matches the EN 1993-1-2 step method', () => {
  // 10 mm steel plate heated on both faces by ISO 834: Am/V = 2/t = 200 1/m; model half thickness with a symmetry plane
  const steel = { thickness: 0.005, k: steelConductivity, rho: 7850, cp: steelSpecificHeat, n: 5 };
  const r = conduction1d({
    layers: [steel], T0C: 20, right: { type: 'convection', h: 25, TinfC: isoFireC, emissivity: 0.7 }, tEnd: 1800, dt: 5,
  });
  const ec3 = steelTemperatureUnprotected({ sectionFactor: 200, fire: 'iso834', duration: 1800, dt: 1 });
  near(r.finalC[0], ec3.finalC, 6, 'FV solution vs lumped step equation at 30 min (273 vs 273.15 offset, implicit vs explicit)');
  assert.ok(r.finalC[5] - r.finalC[0] < 2, 'thin plate is nearly isothermal');
  assert.ok(r.finalC[5] < isoFireC(1800));
});
