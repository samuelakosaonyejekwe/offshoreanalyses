// Indirect evidence on the seabed surface where no survey map or sample covers the site: which grains the waves
// and currents of the site can move (hydrodynamic sorting), how that compares with a mapped sediment class, and
// the confidence grade of the final statement. Pure functions; the relations come from the coastal engine.
import { waveOrbital, bedShearWave, bedShearCurrent, combinedShearSoulsby, criticalShearStress } from '../engine/coastal/sediment.js';

const NU = 1.36e-6; // sea water at about 10 °C, m²/s
/** Finest grain the threshold curve is used for: below very fine sand, cohesion governs and the curve does not apply. */
export const FINEST_SAND = 6.3e-5;
const COARSEST = 0.064;
/** Mean speed over the depth as a fraction of the surface speed for a 1/7 power-law profile. */
const DEPTH_MEAN = 7 / 8;

/**
 * Largest bed shear stress in a wave cycle on a flat bed of grain size d (skin friction, roughness z0 = d/12):
 * linear-theory orbital velocity, Soulsby (1997) rough-bed wave friction factor capped at 0.3, logarithmic
 * current profile, and the Soulsby (1995) wave–current combination with waves and current in line.
 * @param {{ depth: number, Hs?: number, Tp?: number, U?: number, d: number }} o depth [m]; Hs [m]; Tp [s];
 *   U [m/s] depth-mean current; d [m]
 * @returns {{ Uw: number, tauW: number, tauC: number, tauMax: number }} [m/s], [Pa]
 */
export function bedStress({ depth, Hs = 0, Tp = 0, U = 0, d }) {
  const z0 = d / 12;
  let Uw = 0, tauW = 0;
  if (Hs > 0 && Tp > 0) {
    Uw = waveOrbital({ H: Hs, T: Tp, h: depth }).Uw;
    const w = bedShearWave({ Uw, T: Tp, z0 });
    tauW = w.fw > 0.3 ? 0.5 * 1025 * 0.3 * Uw * Uw : w.tau;
  }
  const tauC = U > 0 && depth / z0 > 3 ? bedShearCurrent({ U, h: depth, z0 }).tau : 0;
  return { Uw, tauW, tauC, tauMax: combinedShearSoulsby({ tauC, tauW }).tauMax };
}

/**
 * Largest grain the given conditions set in motion: the coarsest size, between very fine sand and 64 mm, whose
 * Soulsby–Whitehouse threshold stress is reached on a bed of that grain. Null when no sand or gravel size moves.
 * (A rougher bed draws more stress from the same flow, so the sizes are scanned from the coarse end.)
 * @returns {number|null} grain size [m]
 */
export function largestMobileGrain(cond) {
  const mobile = (d) => bedStress({ ...cond, d }).tauMax >= criticalShearStress({ d, nu: NU }).tauCr;
  const n = 60, at = (i) => FINEST_SAND * (COARSEST / FINEST_SAND) ** (i / n);
  for (let i = n; i >= 0; i--) {
    if (!mobile(at(i))) continue;
    if (i === n) return COARSEST;
    let lo = Math.log(at(i)), hi = Math.log(at(i + 1)); // mobile at lo, still at hi
    for (let k = 0; k < 30; k++) { const mid = (lo + hi) / 2; if (mobile(Math.exp(mid))) lo = mid; else hi = mid; }
    return Math.exp(lo);
  }
  return null;
}

const mm = (d) => `${Number((d * 1e3).toPrecision(2))} mm`;
const sizeName = (d) => (d >= 0.002 ? 'gravel' : d >= 5e-4 ? 'coarse sand' : d >= 2.5e-4 ? 'medium sand' : d >= 1.25e-4 ? 'fine sand' : 'very fine sand');

/** Fine sand: the grain size whose movement in frequent conditions marks a bed that is worked all the time. */
export const FINE_SAND = 2e-4;
/**
 * Bed shear stress [Pa] below which mud settles and stays: the critical stress for deposition of cohesive
 * sediment is of the order of 0.05–0.1 Pa (Krone 1962; Whitehouse, Soulsby, Roberts & Mitchener 2000, Dynamics
 * of Estuarine Muds); the lower end is used.
 */
export const MUD_DEPOSITION_STRESS = 0.05;

/**
 * Hydrodynamic sorting at the site. Two states are examined: frequent conditions (the wave height and the
 * current speed that are each exceeded 10 % of the time) and the 1-year storm (1-year wave height and period with
 * the largest current of the record). The outcome is one of three settings:
 *   'worked'  — the stress of frequent conditions on a fine-sand bed reaches the threshold of motion of fine sand (0.2 mm): fine material cannot settle for long, so the bed is
 *               sand or coarser unless a river supplies mud faster than the sea removes it;
 *   'calm'    — in frequent conditions the stress on the bed stays below the level at which mud settles and
 *               stays: mud accumulates wherever fine sediment is supplied, and sand found here is a relict;
 *   'between' — neither: sand, muddy sand and mud all occur in such settings.
 * Peak periods follow Tp = tpFactor·√Hs, the relation the site pull fits to its storms (4.4 when it has none).
 * @param {{ depth: number, waves?: { p90: number, tpFactor?: number, hs1: number, tp1?: number },
 *   current?: { p90: number, max: number } }} o depth [m]; wave heights [m]; periods [s]; surface current [m/s]
 * @returns {null|{ setting: string, frequent: object, storm: object, dFrequent: number|null, dStorm: number|null, text: string }}
 */
export function hydroSorting({ depth, waves, current }) {
  if (!(depth > 0) || (!(waves?.p90 > 0) && !(current?.p90 > 0))) return null;
  const tpOf = (hs, tp) => (tp > 0 ? tp : (waves?.tpFactor || 4.4) * Math.sqrt(hs));
  const frequent = { depth, Hs: waves?.p90 || 0, Tp: waves?.p90 > 0 ? tpOf(waves.p90) : 0, U: DEPTH_MEAN * (current?.p90 || 0) };
  const storm = { depth, Hs: waves?.hs1 || 0, Tp: waves?.hs1 > 0 ? tpOf(waves.hs1, waves.tp1) : 0, U: DEPTH_MEAN * (current?.max || 0) };
  const dFrequent = largestMobileGrain(frequent), dStorm = largestMobileGrain(storm);
  // stresses are quoted on a fine-sand bed so that sites can be compared
  const st = (c) => { const s = bedStress({ ...c, d: FINE_SAND }); return { ...c, Uw: s.Uw, tau: s.tauMax }; };
  const f = st(frequent), s = st(storm);
  const tauFine = criticalShearStress({ d: FINE_SAND, nu: NU }).tauCr;
  const setting = f.tau >= tauFine ? 'worked' : f.tau < MUD_DEPOSITION_STRESS ? 'calm' : 'between';
  const stormText = dStorm ? `the 1-year storm moves grains up to ${mm(dStorm)}` : 'not even the 1-year storm moves very fine sand';
  const text = setting === 'worked'
    ? `Waves and currents that occur about one day in ten put ${f.tau.toFixed(2)} Pa on the bed, enough to move fine sand (${tauFine.toFixed(2)} Pa)${dFrequent ? ` and grains up to ${mm(dFrequent)} (${sizeName(dFrequent)})` : ''}, and ${stormText}: the bed is worked all the time, so mud does not stay unless a river supplies it faster than the sea removes it.`
    : setting === 'calm'
      ? `Even in the conditions of the roughest day in ten, the stress on the bed (${f.tau.toFixed(3)} Pa) stays below the ${MUD_DEPOSITION_STRESS} Pa at which mud settles and stays (${stormText}): fine sediment accumulates wherever it is supplied, and any sand is a relict of an earlier sea level.`
      : `In the conditions of the roughest day in ten the stress on the bed is ${f.tau.toFixed(2)} Pa — too much for mud to settle undisturbed, too little to move fine sand (${tauFine.toFixed(2)} Pa) — and ${stormText}: sand, muddy sand and mud are all found in such settings.`;
  return { setting, frequent: f, storm: s, dFrequent, dStorm, text };
}

// what each mapped class is, for the comparison with the hydrodynamic setting
const FAMILY = {
  rock: 'rock', gravel: 'sandy', sand: 'sandy', muddySand: 'sandy', biogenic: 'sandy', mixed: 'mixed', volcanic: 'mixed',
  sandyMud: 'mud', mud: 'mud', silt: 'mud', clay: 'mud', calcOoze: 'mud', silOoze: 'mud', mixedOoze: 'mud', silMud: 'mud', calcMud: 'mud', lithogenous: 'mud',
};
const SAMPLE = 'Look for a bottom note on the nautical chart or a published survey of the area, or take a grab sample.';

/**
 * Compare a mapped sediment class with the hydrodynamic setting.
 *   'agree'    — two independent lines point the same way: the measured grain size of the class may be passed on;
 *   'conflict' — they contradict each other: nothing is passed on, both are shown with the likely explanations;
 *   'open'     — the hydrodynamics neither support nor contradict the class (or could not be worked out).
 * @param {{ cls: string, d50?: number|null, d50Range?: number[]|null, hydro: object|null }} o measured typical
 *   grain size of the class and its range [m]
 * @returns {{ verdict: string, d50: number|null, range: number[]|null, why: string, check: string[] }}
 */
export function combineSeabedEvidence({ cls, d50 = null, d50Range = null, hydro }) {
  const fam = FAMILY[cls], open = (why) => ({ verdict: 'open', d50: null, range: null, why, check: [] });
  if (!hydro) return open('The wave and current data needed for the hydrodynamic check were not available.');
  if (!fam || fam === 'rock') return open('The hydrodynamic check says nothing about rock.');
  if (fam === 'mixed') return open('A mixed bed can take any grain size, so the hydrodynamics cannot confirm one.');
  if (hydro.setting === 'between') return open('Sand, muddy sand and mud all occur where the bed is neither worked all the time nor left at rest, so the hydrodynamics neither confirm nor contradict the map.');
  if (hydro.setting === 'worked') {
    if (fam === 'sandy') return { verdict: 'agree', d50: d50 == null ? null : Math.max(d50, FINEST_SAND), range: d50Range && [Math.max(d50Range[0], FINEST_SAND), Math.max(d50Range[1], FINEST_SAND)],
      why: 'The map shows a sandy or coarser bed, and frequent waves and currents are strong enough to keep mud from settling.', check: [] };
    return { verdict: 'conflict', d50: null, range: null, why: 'The map shows fine sediment, yet frequent waves and currents here are strong enough to move fine sand.',
      check: ['A river-fed mud belt: near a large river mouth mud is supplied faster than the sea removes it, and the map is right.', 'A swept bed: the map is too coarse here and the real seabed is sand, or a stiff older clay laid bare by erosion.', SAMPLE] };
  }
  if (fam === 'mud') return { verdict: 'agree', d50, range: d50Range, why: 'The map shows fine sediment, and the bed is so rarely disturbed that fine sediment stays where it settles.', check: [] };
  return { verdict: 'conflict', d50: null, range: null, why: 'The map shows sand or coarser sediment, yet the bed here is so rarely disturbed that mud would settle and stay.',
    check: ['Relict sand or gravel: left by a lower sea level or by ice and never buried — common on outer shelves, and then the map is right.', 'Sand made on the spot from shell and coral (carbonate shelves).', 'A thin cover of mud over the mapped sediment, which a coarse map does not show.', SAMPLE] };
}

/** Confidence grades of the seabed statement, best first. */
export const CONFIDENCE = {
  observed: { label: 'Your own observation at the site', chip: 'ok', raise: 'A grab sample describes the top few centimetres at one spot: a geotechnical investigation (boreholes, cone penetration tests) is still needed for strength and layering.' },
  measured: { label: 'Measured samples nearby', chip: 'ok', raise: 'Samples taken at the site itself, and a geotechnical investigation for what lies below the surface.' },
  mapped: { label: 'Mapped by a geological survey', chip: 'ok', raise: 'A sample or survey at the site: a map polygon is drawn from samples that may be kilometres apart.' },
  consistent: { label: 'Consistent indirect evidence', chip: 'info', raise: 'A nautical chart bottom note, a published survey of the area or a single grab sample would turn this into site knowledge.' },
  single: { label: 'Single coarse source', chip: 'warn', raise: 'Any second line of evidence: a chart bottom note, a published survey nearby, or a grab sample.' },
  conflict: { label: 'Sources disagree', chip: 'fail', raise: 'A grab sample or a published survey of the area: until then treat the seabed type as unknown.' },
  none: { label: 'No information', chip: 'fail', raise: 'Any description of the seabed: a chart bottom note, a regional geological map or a sample.' },
};
/** Grade of a `data.seabed` record. */
export function seabedConfidence(sb) {
  const key = !sb || !sb.cls ? 'none' : sb.basis === 'observation' ? 'observed' : sb.basis === 'samples' || sb.basis === 'grainsize' ? 'measured' : sb.basis === 'map' ? 'mapped'
    : sb.evidence?.verdict === 'agree' ? 'consistent' : sb.evidence?.verdict === 'conflict' ? 'conflict' : 'single';
  return { key, ...CONFIDENCE[key] };
}
