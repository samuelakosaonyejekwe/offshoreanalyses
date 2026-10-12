// Shared lateral pile model of the suite: translates the soil-reaction inputs (API-type curves from the soil
// profile, user-supplied p–y curves such as PISA-calibrated reactions, or a linear Winkler modulus) into the
// options of the beam-on-nonlinear-Winkler-foundation solver.
import { pyUser } from '../../engine/geotech/py.js';
import { need, tube } from './_common.js';

/** Soil-reaction model inputs shared by the pile modules. */
export function pyInputs({ group = 'Soil reaction model' } = {}) {
  return [
    { key: 'pyMode', label: 'Lateral soil reactions', type: 'select', default: 'api', group, role: 'soil',
      options: [{ value: 'api', label: 'API-type p–y curves from the soil profile' }, { value: 'user', label: 'User p–y curves (e.g. PISA-calibrated)' }, { value: 'linear', label: 'Linear Winkler modulus (constant)' }],
      help: 'API curves were derived for slender piles and misjudge the stiffness of large-diameter monopiles: supply site-calibrated curves where available.' },
    { key: 'loading', label: 'Loading type for the API curves', type: 'select', default: 'static', group, role: 'soil', showIf: { key: 'pyMode', in: ['api'] },
      options: [{ value: 'static', label: 'Static (monotonic)' }, { value: 'cyclic', label: 'Cyclic (degraded)' }],
      help: 'Cyclic curves carry the API degradation of the ultimate resistance near the surface (A = 0.9 in sand; reduced plateau in clay).' },
    { key: 'pMult', label: 'p-multiplier', unit: '-', type: 'number', default: 1, min: 0.1, max: 5, group, role: 'soil', showIf: { key: 'pyMode', in: ['api'] }, uncertain: { dist: 'lognormal', cov: 0.2 },
      help: 'Scales soil resistance: < 1 for cyclic degradation or installation disturbance, > 1 to match a calibrated large-diameter response.' },
    { key: 'yMult', label: 'y-multiplier', unit: '-', type: 'number', default: 1, min: 0.1, max: 10, group, role: 'soil', showIf: { key: 'pyMode', in: ['api'] },
      help: 'Stretches the displacement axis: > 1 softens the response without changing the ultimate resistance.' },
    { key: 'kLin', label: 'Winkler modulus', unit: 'MPa', type: 'number', default: 50, min: 0.1, max: 5000, group, role: 'soil', showIf: { key: 'pyMode', in: ['linear'] },
      help: 'Soil reaction per unit pile length per unit deflection (N/m² = Pa), constant with depth.' },
    { key: 'pyTable', label: 'User p–y points', type: 'table', group, role: 'soil', showIf: { key: 'pyMode', in: ['user'] },
      columns: [{ key: 'z', label: 'Depth below mudline', unit: 'm', type: 'number' }, { key: 'y', label: 'Deflection y', unit: 'mm', type: 'number' }, { key: 'p', label: 'Reaction p', unit: 'kN/m', type: 'number' }],
      default: [
        { z: 0, y: 0, p: 0 }, { z: 0, y: 20, p: 150 }, { z: 0, y: 100, p: 300 }, { z: 0, y: 600, p: 400 },
        { z: 10, y: 0, p: 0 }, { z: 10, y: 20, p: 2500 }, { z: 10, y: 100, p: 6000 }, { z: 10, y: 600, p: 9000 },
        { z: 40, y: 0, p: 0 }, { z: 40, y: 20, p: 9000 }, { z: 40, y: 100, p: 24000 }, { z: 40, y: 600, p: 40000 },
      ],
      help: 'One curve per depth, each starting at y = 0, p = 0; reactions are interpolated linearly between depths and held constant beyond the first and last depth.' },
  ];
}

function userCurves(rows) {
  const by = new Map();
  for (const r of rows || []) {
    const z = Number(r.z), y = Number(r.y) / 1e3, p = Number(r.p) * 1e3;
    if (!Number.isFinite(z) || !Number.isFinite(y) || !Number.isFinite(p)) continue;
    if (!by.has(z)) by.set(z, []);
    by.get(z).push([y, p]);
  }
  const zs = [...by.keys()].sort((a, b) => a - b);
  need(zs.length >= 1, 'The user p–y table is empty: enter at least one curve (three or more points at one depth).');
  const curves = zs.map((z) => {
    const pts = by.get(z).sort((a, b) => a[0] - b[0]);
    need(pts.length >= 2 && pts[0][0] === 0, `User p–y curve at ${z} m must start at y = 0 and have at least two points.`);
    for (let i = 1; i < pts.length; i++) need(pts[i][0] > pts[i - 1][0], `User p–y curve at ${z} m has repeated deflection values.`);
    return pyUser({ y: pts.map((q) => q[0]), p: pts.map((q) => q[1]) });
  });
  return (z) => {
    if (z <= zs[0]) return curves[0];
    if (z >= zs[zs.length - 1]) return curves[curves.length - 1];
    let i = 0;
    while (zs[i + 1] < z) i++;
    const w = (z - zs[i]) / (zs[i + 1] - zs[i]), a = curves[i], b = curves[i + 1];
    return { type: 'user', pu: (1 - w) * a.pu + w * b.pu, p: (y) => (1 - w) * a.p(y) + w * b.p(y), k: (y) => (1 - w) * a.k(y) + w * b.k(y), meta: {} };
  };
}

/**
 * solveLateralPile options (without head loads) from the module inputs.
 * @param {object} inp module inputs (D m, t mm, Lp m, nElem, pyMode, loading, pMult, yMult, kLin MPa, pyTable)
 * @param {object} profile engine soil profile (used by the API mode)
 * @param {object} steel { E }
 * @param {object} [over] overrides (L, loading, axialLoad, …)
 */
export function pileOptions(inp, profile, steel, over = {}) {
  const sec = tube(inp.D, inp.t / 1e3);
  const base = { L: inp.Lp, D: inp.D, t: inp.t / 1e3, E: steel.E, nElem: inp.nElem ?? 40, ...over };
  need(base.L >= 1.5 * inp.D, `Embedded length ${base.L} m is less than 1.5 diameters: the beam-on-springs model is not meaningful for such a stubby pile.`);
  if (inp.pyMode === 'linear') return { ...base, pyAt: () => inp.kLin * 1e6, sec };
  if (inp.pyMode === 'user') return { ...base, pyAt: userCurves(inp.pyTable), sec };
  need(profile.depth >= base.L, `The soil profile (${profile.depth} m) is shallower than the pile toe (${base.L} m): extend the last layer.`);
  return { ...base, profile, loading: over.loading ?? inp.loading ?? 'static', pMult: inp.pMult ?? 1, yMult: inp.yMult ?? 1, sec };
}
