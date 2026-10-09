// Temporary shell probe: cantilever under tip load by finite differences (replaced by the real suite).
export default {
  id: 's01-probe', suite: 1, title: 'Shell probe — cantilever', summary: 'Cantilever tip deflection, used to exercise the application shell.',
  covers: [2], fidelity: 'numerical', standards: ['Euler–Bernoulli beam theory'], models: ['Euler–Bernoulli beam theory'], limits: 'Probe only.',
  inputs: [
    { key: 'L', label: 'Length', unit: 'm', type: 'number', default: 10, min: 0.1, max: 200, group: 'Geometry', role: 'geometry', design: true },
    { key: 'D', label: 'Outer diameter', unit: 'mm', type: 'number', default: 600, min: 50, max: 5000, group: 'Geometry', role: 'geometry', design: true },
    { key: 't', label: 'Wall thickness', unit: 'mm', type: 'number', default: 20, min: 3, max: 150, group: 'Geometry', role: 'geometry', design: true, uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'mat', label: 'Steel', type: 'select', source: 'materials:steel', default: 's355', group: 'Material', role: 'material' },
    { key: 'P', label: 'Tip load', unit: 'kN', type: 'number', default: 100, min: 0, max: 1e5, group: 'Loading', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.15 } },
    { key: 'n', label: 'Elements', type: 'int', default: 8, min: 2, max: 2000, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'n', label: 'Elements', levels: [4, 8, 16, 32, 64], qoi: 'tip' },
  compute(inp, ctx) {
    const m = ctx.material(inp.mat);
    const D = inp.D / 1e3, t = inp.t / 1e3, I = Math.PI / 64 * (D ** 4 - (D - 2 * t) ** 4), E = m.E, P = inp.P * 1e3, L = inp.L;
    // trapezoidal double integration of M/EI: second-order accurate in 1/n
    const n = inp.n, h = L / n; let th = 0, y = 0; const xs = [0], ys = [0];
    for (let i = 0; i < n; i++) { const k0 = P * (L - i * h) / (E * I), k1 = P * (L - (i + 1) * h) / (E * I); const th1 = th + 0.5 * (k0 + k1) * h; y += 0.5 * (th + th1) * h; th = th1; xs.push((i + 1) * h); ys.push(y * 1e3); }
    const sig = P * L * (D / 2) / I;
    return {
      kpis: [{ key: 'tip', label: 'Tip deflection', value: y * 1e3, unit: 'mm' }, { key: 'sig', label: 'Root bending stress', value: sig / 1e6, unit: 'MPa', limit: m.fy / 1e6 / 1.15, sense: 'max' }],
      checks: [{ name: 'Root bending', demand: sig / 1e6, capacity: m.fy / 1e6 / 1.15, unit: 'MPa', uc: sig / (m.fy / 1.15), ref: 'fy/γm' }],
      plots: [{ type: 'line', title: 'Deflected shape', x: { label: 'x', unit: 'm' }, y: { label: 'Deflection', unit: 'mm' }, series: [{ name: 'w', x: xs, y: ys }] }],
      notes: ['Probe module.'], publish: { 'jacket.baseShear': P }, quantities: { steel_t: Math.PI * (D - t) * t * L * 7.85 },
    };
  },
  verify: [{ name: 'PL³/3EI', source: 'closed form', inputs: { n: 400 }, expect: [{ kpi: 'tip', value: 100e3 * 1e3 / (3 * 210e9 * Math.PI / 64 * (0.6 ** 4 - 0.56 ** 4)) * 1e3, tol: 0.001 }] }],
};
