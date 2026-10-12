// Extreme environmental loading: wind, current and mean wave-drift force on the floater for the design storm
// swept over direction, the resulting offsets and line tensions, and a return-period extrapolation.
import { dist } from '../../engine/math/stats.js';
import { DEG, need, round, fmt, adv, check, worst, sev, systemInputs, buildSystem, fastSystem, lineTable, steelAdvice, LIMIT_QS, LIMIT_PROPS, CONFIRM } from './_common.js';
import { envInputs, hullInputs, envForce, LIMIT_ENV } from './_vessel.js';

export default {
  id: 's10-extreme-environment',
  suite: 10,
  title: 'Extreme environmental loading',
  summary: 'Builds the mean wind, current and wave-drift force for the design storm from every direction and finds the direction and return period that govern offset and line tension.',
  covers: [5],
  fidelity: 'screening',
  standards: ['API RP 2SK (3rd ed., 2005) §4 (environmental criteria) and §5.3', 'DNV-OS-E301 (2018) Ch.2 Sec.1 (environmental conditions and loads)', 'DNV-RP-C205 (2019) §5–7', 'ISO 19901-1 (2015) metocean design criteria'],
  models: ['Wind and current drag on projected areas (cross-flow principle)', 'JONSWAP wave spectrum', 'Mean wave-drift force 2∫S·D dω with a Maruo reflection coefficient', 'Quasi-static mooring equilibrium per direction', 'Gumbel extreme-value extrapolation over return period'],
  limits: `${LIMIT_ENV} ${LIMIT_QS} Only mean loads are applied here: low- and wave-frequency dynamics are added in the tension-and-offset module. The return-period extrapolation fits a Gumbel law through two return levels and is indicative beyond them. ${LIMIT_PROPS}`,
  publishes: ['vessel.Fmean'],
  searchTerms: 'mooring extreme environment 100-year wind wave current directional',
  inputs: [
    ...systemInputs(),
    ...hullInputs(),
    ...envInputs(),
    { key: 'Hs1', label: '1-year significant wave height', unit: 'm', type: 'number', default: 7.5, min: 0, max: 20, group: 'Environment (1-year, for extrapolation)', role: 'environment', bind: 'site.Hs1' },
    { key: 'Uw1', label: '1-year wind speed', unit: 'm/s', type: 'number', default: 25, min: 0, max: 60, group: 'Environment (1-year, for extrapolation)', role: 'environment', bind: 'site.Uw1' },
    { key: 'Uc1', label: '1-year current speed', unit: 'm/s', type: 'number', default: 0.7, min: 0, max: 3, group: 'Environment (1-year, for extrapolation)', role: 'environment', bind: 'site.Uc1' },
    { key: 'returnYears', label: 'Return period of the design storm', unit: 'years', type: 'number', default: 100, min: 5, max: 10000, group: 'Environment (design storm)', role: 'criteria',
      help: 'Return period attached to the design-storm values above (100 years for permanent moorings).' },
    { key: 'waveOffDeg', label: 'Wave direction relative to wind', unit: 'deg', type: 'number', default: 0, min: -90, max: 90, group: 'Directional combination', role: 'bc', help: 'Collinear (0) is the usual governing case for spread moorings.' },
    { key: 'curOffDeg', label: 'Current direction relative to wind', unit: 'deg', type: 'number', default: 0, min: -180, max: 180, group: 'Directional combination', role: 'bc', help: 'API RP 2SK also asks for a crossed case (e.g. current at 30–45° to wind and waves) for weathervaning units.' },
    { key: 'headingDeg0', label: 'Vessel heading', unit: 'deg', type: 'number', default: 0, min: -180, max: 180, group: 'Directional combination', role: 'ic', help: 'Heading of the hull x-axis, held fixed during the sweep (initial orientation of the vessel).' },
    { key: 'nDir', label: 'Directions in the sweep', type: 'int', default: 16, min: 4, max: 72, group: 'Directional combination', role: 'control' },
    { key: 'sfStatic', label: 'Required safety factor on mean tension', unit: '-', type: 'number', default: 2.0, min: 1, max: 5, group: 'Acceptance criteria', role: 'criteria', help: `Quasi-static screen on the mean tension (API RP 2SK intact, quasi-static). ${CONFIRM}` },
    { key: 'offsetLimitPct', label: 'Allowable mean offset', unit: '% depth', type: 'number', default: 5, min: 0.5, max: 60, group: 'Acceptance criteria', role: 'criteria' },
  ],
  compute(inp) {
    const b = buildSystem(inp);
    const fs = fastSystem(b.system, { maxOffset: Math.max(0.12 * b.depth, 25) });
    const calm = fs.solve([0, 0, 0, 0, 0, 0]);
    const at = (env, dirDeg) => {
      const f = envForce({ ...inp, ...env }, { windDir: dirDeg, waveDir: dirDeg + inp.waveOffDeg, curDir: dirDeg + inp.curOffDeg, headingDeg: inp.headingDeg0, rhoW: inp.rhoW });
      const eq = fs.solve([f.Fx, f.Fy, 0, 0, 0, 0], [calm.x, calm.y]);
      return { f, eq };
    };
    const design = { Hs: inp.Hs, Uw: inp.Uw, Uc: inp.Uc };
    const rows = [];
    for (let k = 0; k < inp.nDir; k++) {
      const dir = 360 * k / inp.nDir;
      const { f, eq } = at(design, dir);
      rows.push({ dir, f, eq });
    }
    const gT = rows.reduce((a, r) => (r.eq.maxTension > a.eq.maxTension ? r : a), rows[0]);
    const gO = rows.reduce((a, r) => (r.eq.offset > a.eq.offset ? r : a), rows[0]);
    const gF = rows.reduce((a, r) => (r.f.F > a.f.F ? r : a), rows[0]);
    const head = envForce(inp, { windDir: 0, waveDir: 0, curDir: 0, headingDeg: 0, rhoW: inp.rhoW });
    // return-period extrapolation of the governing-direction tension: Gumbel through the 1-year value (taken as the mode of
    // the annual maximum, non-exceedance e^-1) and the design-storm value at its return period
    const one = at({ Hs: Math.min(inp.Hs1, inp.Hs), Uw: Math.min(inp.Uw1, inp.Uw), Uc: Math.min(inp.Uc1, inp.Uc) }, gT.dir);
    const yT = -Math.log(-Math.log(1 - 1 / inp.returnYears));
    const gum = (v1, vT) => { const scale = Math.max((vT - v1) / yT, 1e-9 * Math.abs(vT) + 1e-12); return (T) => dist({ dist: 'gumbel', loc: v1, scale }).icdf(1 - 1 / T); };
    const Tof = gum(one.eq.maxTension, gT.eq.maxTension), Fof = gum(one.f.F, gT.f.F);
    const periods = [5, 10, 50, 100, 1000, 10000];
    const util = gT.eq.maxTension / b.MBL;
    const checks = [
      check('Mean tension in the governing direction', gT.eq.maxTension * inp.sfStatic / 1e3, b.MBL / 1e3, 'kN', 'API RP 2SK (3rd ed.) §6.3 (quasi-static factor on mean load, screening)',
        'Guards against line breakage under the steady part of the design storm. Dynamic tension is not included, so this is a necessary screen before the full tension check.'),
      check('Mean offset in the governing direction', gO.eq.offset, inp.offsetLimitPct / 100 * b.depth, 'm', 'API RP 2SK (3rd ed.) §6.2; project riser limit',
        'Guards against riser or umbilical damage from the steady excursion in the worst direction.'),
      { name: 'Anchor uplift in the governing direction', demand: round(Math.max(...gT.eq.lines.map((l) => l.anchorV)) / 1e3, 4), capacity: 0, unit: 'kN', uc: Math.max(...gT.eq.lines.map((l) => l.anchorV)) > 1 ? 99 : 0, ref: 'API RP 2SK (3rd ed.) §6.5',
        explain: 'A non-zero vertical anchor load means the windward line has lifted completely off the seabed under the mean storm load; drag anchors cannot take it.' },
    ];
    const g = worst(checks);
    const T10k = Tof(10000);
    return {
      kpis: [
        { key: 'Fenv', label: 'Largest mean environmental force', value: round(gF.f.F / 1e3, 5), unit: 'kN', explain: 'Resultant of wind, current and mean wave drift in the direction where it is largest. This is the mean load the mooring must balance.' },
        { key: 'dirF', label: 'Direction of largest force', value: round(gF.dir, 4), unit: 'deg', explain: 'Direction the environment travels towards when the force is largest (beam-on for elongated hulls).' },
        { key: 'FwindHead', label: 'Wind force, head-on', value: round(head.wind / 1e3, 5), unit: 'kN', explain: 'Wind drag ½·ρair·Cd·A·U² on the frontal projected area; usually the largest steady load on a semi-submersible.' },
        { key: 'FcurHead', label: 'Current force, head-on', value: round(head.current / 1e3, 5), unit: 'kN', explain: 'Current drag ½·ρ·Cd·A·U² on the frontal underwater area.' },
        { key: 'FdriftHead', label: 'Mean wave-drift force, head-on', value: round(head.drift / 1e3, 5), unit: 'kN', explain: 'Steady second-order wave force from reflection of wave energy by the hull, integrated over the spectrum.' },
        { key: 'offsetMax', label: 'Largest mean offset', value: round(gO.eq.offset, 5), unit: 'm', limit: round(inp.offsetLimitPct / 100 * b.depth, 5), sense: 'max', explain: 'Mean offset in the worst direction (normally a load between two lines, where the system is softest).' },
        { key: 'dirOffset', label: 'Direction of largest offset', value: round(gO.dir, 4), unit: 'deg', explain: 'Environmental direction giving the largest mean offset.' },
        { key: 'Tmax', label: 'Largest mean line tension', value: round(gT.eq.maxTension / 1e3, 5), unit: 'kN', explain: 'Mean tension of the most loaded line in the worst direction (normally a load in line with one line).' },
        { key: 'dirT', label: 'Direction of largest tension', value: round(gT.dir, 4), unit: 'deg', explain: 'Environmental direction giving the largest line tension; use it as the intact governing case in the tension-and-offset module.' },
        { key: 'util', label: 'Largest mean tension / MBL', value: round(util, 4), unit: '-', limit: round(1 / inp.sfStatic, 4), sense: 'max', explain: 'Share of the break load used by the mean storm tension.' },
        { key: 'T10000', label: 'Mean tension extrapolated to 10 000 years', value: round(T10k / 1e3, 5), unit: 'kN', explain: 'Gumbel extrapolation of the governing mean tension to the 10 000-year level, an indication of robustness (abnormal environment) rather than a design value.' },
        { key: 'ratio1', label: '1-year / design-storm tension', value: round(one.eq.maxTension / gT.eq.maxTension, 4), unit: '-', explain: 'How much of the design-storm mean tension is already present in the annual storm; a high ratio means pretension dominates and fatigue will matter.' },
      ],
      checks,
      plots: [
        { type: 'polar', title: 'Mean environmental force by direction', series: [{ name: 'Total (kN)', thetaDeg: rows.map((r) => r.dir), r: rows.map((r) => round(r.f.F / 1e3)) }, { name: 'Wind (kN)', thetaDeg: rows.map((r) => r.dir), r: rows.map((r) => round(r.f.wind / 1e3)) },
          { name: 'Current (kN)', thetaDeg: rows.map((r) => r.dir), r: rows.map((r) => round(r.f.current / 1e3)) }, { name: 'Wave drift (kN)', thetaDeg: rows.map((r) => r.dir), r: rows.map((r) => round(r.f.drift / 1e3)) }],
          explain: 'Force rose: how the steady load and its three parts vary with the direction the storm travels towards. Lobes show the beam-on directions where the hull presents most area.' },
        { type: 'polar', title: 'Mean offset and largest line tension by direction', series: [{ name: 'Offset (m × 100)', thetaDeg: rows.map((r) => r.dir), r: rows.map((r) => round(r.eq.offset * 100)) }, { name: 'Largest tension (kN)', thetaDeg: rows.map((r) => r.dir), r: rows.map((r) => round(r.eq.maxTension / 1e3)) }],
          explain: 'Response rose. Offset peaks between lines and tension peaks in line with a line; with n lines both patterns repeat n times around the compass. Offsets are scaled by 100 to share the axis.' },
        { type: 'line', title: 'Return-period extrapolation (governing direction)', x: { label: 'Return period', unit: 'years (log10)' }, y: { label: 'Force / tension', unit: 'kN' },
          series: [{ name: 'Mean environmental force', x: periods.map((p) => round(Math.log10(p))), y: periods.map((p) => round(Fof(p) / 1e3)) }, { name: 'Largest mean line tension', x: periods.map((p) => round(Math.log10(p))), y: periods.map((p) => round(Tof(p) / 1e3)) },
            { name: 'MBL / safety factor', x: [Math.log10(5), 4], y: [b.MBL / inp.sfStatic / 1e3, b.MBL / inp.sfStatic / 1e3], dash: true }],
          explain: 'Gumbel lines through the 1-year and design-storm values, against log return period. Where the tension line meets the allowable is the return period at which the mean-load screen is used up.' },
      ],
      tables: [
        { title: 'Directional results for the design storm', columns: [{ key: 'dir', label: 'Direction (towards)', unit: 'deg' }, { key: 'F', label: 'Total force', unit: 'kN' }, { key: 'wind', label: 'Wind', unit: 'kN' }, { key: 'cur', label: 'Current', unit: 'kN' }, { key: 'drift', label: 'Wave drift', unit: 'kN' }, { key: 'off', label: 'Mean offset', unit: 'm' }, { key: 'T', label: 'Largest tension', unit: 'kN' }, { key: 'line', label: 'Critical line' }],
          rows: rows.map((r) => ({ dir: round(r.dir), F: round(r.f.F / 1e3), wind: round(r.f.wind / 1e3), cur: round(r.f.current / 1e3), drift: round(r.f.drift / 1e3), off: round(r.eq.offset), T: round(r.eq.maxTension / 1e3), line: b.system.lines[r.eq.maxTensionLine].name })),
          explain: 'One row per storm direction with the load build-up, the mean offset and which line is critical. Use it to pick the intact and damaged governing cases.' },
        { title: 'Return-period levels in the governing direction', columns: [{ key: 'T', label: 'Return period', unit: 'years' }, { key: 'F', label: 'Mean force', unit: 'kN' }, { key: 'Tl', label: 'Largest mean tension', unit: 'kN' }, { key: 'u', label: 'Tension / MBL', unit: '%' }],
          rows: periods.map((p) => ({ T: p, F: round(Fof(p) / 1e3), Tl: round(Tof(p) / 1e3), u: round(100 * Tof(p) / b.MBL, 3) })),
          explain: 'Indicative force and mean tension at other return periods from the two-point Gumbel fit. Values beyond the design return period are extrapolations.' },
        lineTable(gT.eq.lines, b.MBL, 'Line results in the tension-governing direction'),
      ],
      notes: [`Largest force ${fmt(gF.f.F / 1e3)} kN towards ${fmt(gF.dir)}°: wind ${fmt(gF.f.wind / 1e3)}, current ${fmt(gF.f.current / 1e3)}, wave drift ${fmt(gF.f.drift / 1e3)} kN.`,
        `Tension governs towards ${fmt(gT.dir)}° (${b.system.lines[gT.eq.maxTensionLine].name}); offset governs towards ${fmt(gO.dir)}°.`,
        'No current or wave load is applied to the lines themselves; force coefficients are generic inputs.'],
      advice: [
        adv('engineering', sev(g.uc), `The ${inp.returnYears}-year storm gives up to ${fmt(gF.f.F / 1e3)} kN mean force, ${fmt(gO.eq.offset)} m mean offset and ${fmt(gT.eq.maxTension / 1e3)} kN mean tension (${fmt(100 * util)} % MBL). Governing: ${g.name} (${fmt(g.uc)}).`,
          g.uc > 1 ? (g.name.includes('offset') ? `Stiffen the system in the ${fmt(gO.dir)}° sector (more pretension, or rotate the pattern so a line group faces the dominant storm direction).` : g.name.includes('uplift') ? 'Add ground chain or use uplift-capable anchors on the windward lines.' : `Strengthen or add lines facing ${fmt((gT.dir + 180) % 360)}° (the windward side).`)
            : `Carry the ${fmt(gT.dir)}° case (tension) and the ${fmt(gO.dir)}° case (offset) into the tension-and-offset and line-break modules; also run a crossed-current case.`,
          'API RP 2SK (3rd ed., 2005) §4.3 (design environments and directional combinations); DNV-OS-E301 (2018) Ch.2 Sec.1 [3].'),
        adv('engineering', 'info', `Wind contributes ${fmt(100 * gF.f.wind / gF.f.F)} %, current ${fmt(100 * gF.f.current / gF.f.F)} % and wave drift ${fmt(100 * gF.f.drift / gF.f.F)} % of the governing force (scalar shares).`,
          'Spend the calibration effort on the largest contributor: wind-tunnel coefficients for wind, measured current profiles and line drag for current, diffraction-analysis drift coefficients for waves.',
          'DNV-RP-C205 (2019) §5–7; ISO 19901-7 (2013) §7.'),
        steelAdvice(b, Math.max(checks[0].uc, 0.01), 'mean storm-tension'),
      ],
      publish: { 'vessel.Fmean': gF.f.F },
      quantities: { steel_t: round(b.steel_t) },
    };
  },
  verify: [
    { // Head-on wind: F = ½·ρair·Cd·A·U² = 0.5·1.226·1.0·1500·40² = 1471.2 kN.
      name: 'Wind force on the frontal projected area',
      source: 'DNV-RP-C205 (2019) §5.3: F = C·q·S with q = ½·ρ·U²',
      inputs: { Uw: 40, Afront: 1500, CdWind: 1.0 }, expect: [{ kpi: 'FwindHead', value: 1471.2, tol: 0.001 }] },
    { // Deep-draught wall (draught 250 m so all wave energy is intercepted): F = 2∫S·½ρgR²B dω = ρ·g·R²·B·m0 with
      // m0 = Hs²/16 = 6.25 m² for Hs = 10 m → 1025·9.80665·0.25·70·6.25 = 1099.4 kN.
      name: 'Mean wave-drift force of a reflecting wall in irregular waves',
      source: 'Maruo (1960) / Faltinsen (1990) Sea Loads Eq. 5.38 integrated over the spectrum',
      inputs: { Hs: 10, Tp: 12, gammaJ: 3.3, R2: 0.25, Bv: 70, draught: 250 }, expect: [{ kpi: 'FdriftHead', value: 1099.4, tol: 0.02 }] },
    { // Current: F = ½·ρ·Cd·A·U² = 0.5·1025·1.0·1200·1.5² = 1383.8 kN.
      name: 'Current force on the frontal underwater area',
      source: 'DNV-RP-C205 (2019) §6 (drag on projected area)',
      inputs: { Uc: 1.5, AcFront: 1200, CdCur: 1.0 }, expect: [{ kpi: 'FcurHead', value: 1383.8, tol: 0.001 }] },
  ],
};
