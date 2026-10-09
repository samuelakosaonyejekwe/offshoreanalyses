// Cycle counting of load and stress histories. All functions are unit-agnostic: ranges, means and levels come out
// in the unit of the input signal.

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function isSeries(x) {
  return Array.isArray(x) || ArrayBuffer.isView(x);
}

/**
 * Peak–valley (reversal) extraction with a hysteresis (racetrack) gate. A turning point is accepted only once
 * the signal has moved away from it by more than `gate`; the first sample and the last tracked extreme are kept
 * so that the residue of open half cycles is preserved.
 * Reference: ASTM E1049-85 §3 (reversal definition); hysteresis filtering as in common rainflow practice.
 * @param {ArrayLike<number>} x signal samples [signal unit]
 * @param {object} [o]
 * @param {number} [o.gate=0] hysteresis threshold [signal unit]; 0 keeps every strict reversal
 * @returns {{values: number[], indices: number[]}} reversal values [signal unit] and their sample indices [-]
 */
export function extractReversals(x, { gate = 0 } = {}) {
  need(isSeries(x), 'x must be an array');
  need(gate >= 0, 'gate must be non-negative');
  const values = [];
  const indices = [];
  const n = x.length;
  if (n === 0) return { values, indices };
  let ext = x[0];
  let extIdx = 0;
  let dir = 0;
  values.push(ext);
  indices.push(0);
  for (let i = 1; i < n; i++) {
    const v = x[i];
    need(Number.isFinite(v), `x[${i}] must be finite`);
    if (dir === 0) {
      if (Math.abs(v - ext) > gate) {
        dir = v > ext ? 1 : -1;
        ext = v;
        extIdx = i;
      }
    } else if (dir > 0) {
      if (v > ext) {
        ext = v;
        extIdx = i;
      } else if (ext - v > gate) {
        values.push(ext);
        indices.push(extIdx);
        dir = -1;
        ext = v;
        extIdx = i;
      }
    } else if (v < ext) {
      ext = v;
      extIdx = i;
    } else if (v - ext > gate) {
      values.push(ext);
      indices.push(extIdx);
      dir = 1;
      ext = v;
      extIdx = i;
    }
  }
  if (dir !== 0) {
    values.push(ext);
    indices.push(extIdx);
  }
  return { values, indices };
}

function cycle(from, to, count) {
  return { range: Math.abs(to - from), mean: 0.5 * (from + to), count, from, to };
}

/**
 * Rainflow counting by the three-point rule with half cycles.
 * Reference: ASTM E1049-85 (2017) §5.4.4 "Rainflow counting". Ranges that contain the starting point and the
 * ranges left on the stack at the end are counted as half cycles.
 * @param {ArrayLike<number>} x signal samples (reversals are extracted internally) [signal unit]
 * @param {object} [o]
 * @param {number} [o.gate=0] hysteresis gate passed to {@link extractReversals} [signal unit]
 * @returns {{cycles: {range:number, mean:number, count:number, from:number, to:number}[], reversals: number[]}}
 *   count is 1 or 0.5; range, mean, from, to in signal units
 */
export function rainflowThreePoint(x, { gate = 0 } = {}) {
  const rev = extractReversals(x, { gate }).values;
  const cycles = [];
  const stack = [];
  for (let i = 0; i < rev.length; i++) {
    stack.push(rev[i]);
    while (stack.length >= 3) {
      const n = stack.length;
      const xr = Math.abs(stack[n - 1] - stack[n - 2]);
      const yr = Math.abs(stack[n - 2] - stack[n - 3]);
      if (xr < yr) break;
      if (n === 3) {
        cycles.push(cycle(stack[0], stack[1], 0.5));
        stack.shift();
      } else {
        cycles.push(cycle(stack[n - 3], stack[n - 2], 1));
        const last = stack[n - 1];
        stack.length = n - 3;
        stack.push(last);
      }
    }
  }
  for (let i = 0; i + 1 < stack.length; i++) cycles.push(cycle(stack[i], stack[i + 1], 0.5));
  return { cycles, reversals: rev };
}

function fourPointCore(rev, cycles) {
  const stack = [];
  for (let i = 0; i < rev.length; i++) {
    stack.push(rev[i]);
    while (stack.length >= 4) {
      const n = stack.length;
      const inner = Math.abs(stack[n - 2] - stack[n - 3]);
      if (inner > Math.abs(stack[n - 3] - stack[n - 4]) || inner > Math.abs(stack[n - 1] - stack[n - 2])) break;
      cycles.push(cycle(stack[n - 3], stack[n - 2], 1));
      const last = stack[n - 1];
      stack.length = n - 2;
      stack[n - 3] = last;
    }
  }
  return stack;
}

/**
 * Close the cycles left in a four-point residue.
 * 'repeat' (repeat-residue method): the residue is appended to itself and counted again; the cycles that close
 * are exactly those of the history treated as a repeating block. 'half': every residue range is a half cycle
 * (ASTM convention). 'none': nothing is counted.
 * Reference: Amzallag, Gerey, Robert & Bahuaud (1994) Int. J. Fatigue 16, 287–293; Marsh et al. (2016) Int. J.
 * Fatigue 82, 757–765 for the residue treatment.
 * @param {number[]} residue reversal sequence [signal unit]
 * @param {string} [method='repeat'] 'repeat' | 'half' | 'none'
 * @returns {{range:number, mean:number, count:number, from:number, to:number}[]} cycles
 */
export function closeResidue(residue, method = 'repeat') {
  need(isSeries(residue), 'residue must be an array');
  const cycles = [];
  if (method === 'none' || residue.length < 2) return cycles;
  if (method === 'half') {
    for (let i = 0; i + 1 < residue.length; i++) cycles.push(cycle(residue[i], residue[i + 1], 0.5));
    return cycles;
  }
  need(method === 'repeat', `method: unknown residue method "${method}"`);
  const doubled = extractReversals([...residue, ...residue]).values;
  fourPointCore(doubled, cycles);
  return cycles;
}

/**
 * Rainflow counting by the four-point rule (closed hysteresis loops only) with explicit residue.
 * Reference: Amzallag et al. (1994); AFNOR A03-406. With residue = 'half' the result equals the ASTM three-point
 * count; with 'repeat' it equals the ASTM "simplified rainflow for repeating histories" (§5.4.5).
 * @param {ArrayLike<number>} x signal samples [signal unit]
 * @param {object} [o]
 * @param {number} [o.gate=0] hysteresis gate [signal unit]
 * @param {string} [o.residue='repeat'] 'repeat' | 'half' | 'none'
 * @returns {{cycles: object[], residue: number[], reversals: number[]}} cycles as in {@link rainflowThreePoint};
 *   residue is the open reversal sequence before closing
 */
export function rainflowFourPoint(x, { gate = 0, residue = 'repeat' } = {}) {
  const rev = extractReversals(x, { gate }).values;
  const cycles = [];
  const res = fourPointCore(rev, cycles);
  const closing = closeResidue(res, residue);
  for (let i = 0; i < closing.length; i++) cycles.push(closing[i]);
  return { cycles, residue: res, reversals: rev };
}

/**
 * Stitch the residues of consecutive segments of one long history: the residues are concatenated in order,
 * reduced to reversals and counted again, so that the large cycles spanning segment boundaries are recovered.
 * Reference: Marsh et al. (2016) Int. J. Fatigue 82, 757–765.
 * @param {number[][]} residues residue of each segment in chronological order [signal unit]
 * @param {object} [o]
 * @param {string} [o.residue='repeat'] treatment of the final residue: 'repeat' | 'half' | 'none'
 * @returns {{cycles: object[], residue: number[]}} cycles closed by stitching (plus final-residue cycles) and the
 *   final residue before closing
 */
export function stitchResidues(residues, { residue = 'repeat' } = {}) {
  need(Array.isArray(residues), 'residues must be an array of arrays');
  const joined = [];
  for (const r of residues) for (let i = 0; i < r.length; i++) joined.push(r[i]);
  const rev = extractReversals(joined).values;
  const cycles = [];
  const res = fourPointCore(rev, cycles);
  const closing = closeResidue(res, residue);
  for (let i = 0; i < closing.length; i++) cycles.push(closing[i]);
  return { cycles, residue: res };
}

/**
 * Rainflow count of a history supplied in segments (e.g. consecutive measurement files).
 * @param {ArrayLike<number>[]} segments signal segments in chronological order [signal unit]
 * @param {object} [o]
 * @param {number} [o.gate=0] hysteresis gate [signal unit]
 * @param {string} [o.residue='repeat'] treatment of the final residue
 * @returns {{cycles: object[], residue: number[]}}
 */
export function rainflowSegments(segments, { gate = 0, residue = 'repeat' } = {}) {
  need(Array.isArray(segments), 'segments must be an array');
  const cycles = [];
  const residues = [];
  for (const seg of segments) {
    const r = rainflowFourPoint(seg, { gate, residue: 'none' });
    for (let i = 0; i < r.cycles.length; i++) cycles.push(r.cycles[i]);
    residues.push(r.residue);
  }
  const st = stitchResidues(residues, { residue });
  for (let i = 0; i < st.cycles.length; i++) cycles.push(st.cycles[i]);
  return { cycles, residue: st.residue };
}

/**
 * Range histogram of counted cycles.
 * @param {object[]} cycles output of a rainflow function
 * @param {object} [o]
 * @param {number} [o.bins=32] number of equal-width bins [-]
 * @param {number} [o.maxRange] upper edge of the last bin [signal unit]; default = largest range
 * @returns {{edges: Float64Array, centres: Float64Array, counts: Float64Array}} edges (bins+1) and centres in
 *   signal units, counts in cycles
 */
export function rangeHistogram(cycles, { bins = 32, maxRange } = {}) {
  need(Number.isInteger(bins) && bins > 0, 'bins must be a positive integer');
  let top = maxRange;
  if (top === undefined) {
    top = 0;
    for (const c of cycles) if (c.range > top) top = c.range;
  }
  need(top >= 0, 'maxRange must be non-negative');
  if (top === 0) top = 1;
  const w = top / bins;
  const edges = new Float64Array(bins + 1);
  const centres = new Float64Array(bins);
  const counts = new Float64Array(bins);
  for (let i = 0; i <= bins; i++) edges[i] = i * w;
  for (let i = 0; i < bins; i++) centres[i] = (i + 0.5) * w;
  for (const c of cycles) {
    let k = Math.floor(c.range / w);
    if (k >= bins) k = bins - 1;
    counts[k] += c.count;
  }
  return { edges, centres, counts };
}

/**
 * Range–mean matrix of counted cycles.
 * @param {object[]} cycles output of a rainflow function
 * @param {object} [o]
 * @param {number} [o.rangeBins=16] [-]
 * @param {number} [o.meanBins=16] [-]
 * @returns {{rangeEdges: Float64Array, meanEdges: Float64Array, matrix: number[][]}} matrix[iRange][iMean] in
 *   cycles; edges in signal units
 */
export function rangeMeanMatrix(cycles, { rangeBins = 16, meanBins = 16 } = {}) {
  need(Number.isInteger(rangeBins) && rangeBins > 0, 'rangeBins must be a positive integer');
  need(Number.isInteger(meanBins) && meanBins > 0, 'meanBins must be a positive integer');
  let rMax = 0;
  let mMin = Infinity;
  let mMax = -Infinity;
  for (const c of cycles) {
    if (c.range > rMax) rMax = c.range;
    if (c.mean < mMin) mMin = c.mean;
    if (c.mean > mMax) mMax = c.mean;
  }
  if (!(rMax > 0)) rMax = 1;
  if (!(mMax > mMin)) { mMin = (Number.isFinite(mMin) ? mMin : 0) - 0.5; mMax = mMin + 1; }
  const rw = rMax / rangeBins;
  const mw = (mMax - mMin) / meanBins;
  const rangeEdges = Float64Array.from({ length: rangeBins + 1 }, (_, i) => i * rw);
  const meanEdges = Float64Array.from({ length: meanBins + 1 }, (_, i) => mMin + i * mw);
  const matrix = Array.from({ length: rangeBins }, () => new Array(meanBins).fill(0));
  for (const c of cycles) {
    const i = Math.min(rangeBins - 1, Math.floor(c.range / rw));
    const j = Math.min(meanBins - 1, Math.floor((c.mean - mMin) / mw));
    matrix[i][j] += c.count;
  }
  return { rangeEdges, meanEdges, matrix };
}

/**
 * From–to (Markov) matrix of counted cycles: entry [i][j] holds the cycles that start in level class i and
 * turn in level class j.
 * @param {object[]} cycles output of a rainflow function
 * @param {object} [o]
 * @param {number} [o.bins=32] number of level classes [-]
 * @param {number} [o.min] lower bound of the classes [signal unit]; default = smallest turning point
 * @param {number} [o.max] upper bound [signal unit]; default = largest turning point
 * @returns {{edges: Float64Array, matrix: number[][]}}
 */
export function fromToMatrix(cycles, { bins = 32, min, max } = {}) {
  need(Number.isInteger(bins) && bins > 0, 'bins must be a positive integer');
  let lo = min;
  let hi = max;
  if (lo === undefined || hi === undefined) {
    let a = Infinity;
    let b = -Infinity;
    for (const c of cycles) {
      a = Math.min(a, c.from, c.to);
      b = Math.max(b, c.from, c.to);
    }
    if (lo === undefined) lo = Number.isFinite(a) ? a : 0;
    if (hi === undefined) hi = Number.isFinite(b) ? b : 1;
  }
  if (!(hi > lo)) hi = lo + 1;
  const w = (hi - lo) / bins;
  const edges = Float64Array.from({ length: bins + 1 }, (_, i) => lo + i * w);
  const matrix = Array.from({ length: bins }, () => new Array(bins).fill(0));
  const cls = (v) => Math.max(0, Math.min(bins - 1, Math.floor((v - lo) / w)));
  for (const c of cycles) matrix[cls(c.from)][cls(c.to)] += c.count;
  return { edges, matrix };
}

/**
 * Damage-equivalent load: the constant range that, applied nRef times, gives the same Miner damage on a
 * single-slope curve with Wöhler exponent m: DEL = (Σ nᵢ·Sᵢ^m / nRef)^(1/m).
 * Reference: IEC 61400-1 / IEC TS 61400-13 practice; Miner's rule with N ∝ S^(−m).
 * @param {object} o
 * @param {object[]} [o.cycles] counted cycles ({range, count})
 * @param {ArrayLike<number>} [o.ranges] alternative: ranges [signal unit]
 * @param {ArrayLike<number>} [o.counts] counts for `ranges` [cycles]
 * @param {number} o.m Wöhler exponent [-]
 * @param {number} o.nRef reference number of cycles [cycles]
 * @returns {number} damage-equivalent range [signal unit]
 */
export function damageEquivalentLoad({ cycles, ranges, counts, m, nRef }) {
  need(m > 0, 'm must be positive');
  need(nRef > 0, 'nRef must be positive');
  let sum = 0;
  if (cycles) {
    for (const c of cycles) sum += c.count * c.range ** m;
  } else {
    need(isSeries(ranges) && isSeries(counts) && ranges.length === counts.length, 'ranges and counts must be arrays of equal length');
    for (let i = 0; i < ranges.length; i++) sum += counts[i] * ranges[i] ** m;
  }
  return (sum / nRef) ** (1 / m);
}

/**
 * Level-crossing count: number of positive-slope crossings of each level.
 * Reference: ASTM E1049-85 §5.1.
 * @param {ArrayLike<number>} x signal [signal unit]
 * @param {ArrayLike<number>} levels levels [signal unit]
 * @returns {{levels: number[], upCrossings: Float64Array}}
 */
export function levelCrossings(x, levels) {
  need(isSeries(x), 'x must be an array');
  need(isSeries(levels), 'levels must be an array');
  const up = new Float64Array(levels.length);
  for (let k = 0; k < levels.length; k++) {
    const L = levels[k];
    let c = 0;
    for (let i = 1; i < x.length; i++) if (x[i - 1] < L && x[i] >= L) c++;
    up[k] = c;
  }
  return { levels: Array.from(levels), upCrossings: up };
}

/**
 * Total cycle count and largest range of a cycle list.
 * @param {object[]} cycles
 * @returns {{totalCycles:number, maxRange:number}}
 */
export function cycleSummary(cycles) {
  let n = 0;
  let r = 0;
  for (const c of cycles) {
    n += c.count;
    if (c.range > r) r = c.range;
  }
  return { totalCycles: n, maxRange: r };
}
