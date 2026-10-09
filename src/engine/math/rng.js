/**
 * Seeded pseudo-random number generation.
 *
 * Generator: sfc32 (Chris Doty-Humphrey's "Small Fast Counter", PractRand-clean beyond 2^40 outputs), state
 * initialised from the seed through splitmix32 (Steele, Lea & Flood 2014, "Fast splittable pseudorandom number
 * generators", 32-bit finaliser variant) and warmed up. Uniform deviates carry 53 random bits.
 * Normal deviates: Marsaglia polar method (Marsaglia & Bray 1964).
 */

/**
 * splitmix32 stream used to expand a single seed into generator state.
 * @param {number} seed 32-bit integer seed [-]
 * @returns {() => number} function returning successive unsigned 32-bit integers [-]
 */
export function splitmix32(seed) {
  let s = seed | 0;
  return () => {
    s = (s + 0x9e3779b9) | 0;
    let z = s;
    z ^= z >>> 16;
    z = Math.imul(z, 0x21f0aaad);
    z ^= z >>> 15;
    z = Math.imul(z, 0x735a2d97);
    z ^= z >>> 15;
    return z >>> 0;
  };
}

/**
 * Create a deterministic random generator.
 *
 * Reference: sfc32 (PractRand, C. Doty-Humphrey); seeding by splitmix32; normal deviates by the Marsaglia polar
 * method. Validity: any finite seed (non-integers are hashed through their integer and fractional parts);
 * period >= 2^32 guaranteed, ~2^127 on average.
 *
 * @param {number} [seed=1] seed [-]
 * @returns {{ next: () => number, normal: () => number, int: (n: number) => number, shuffle: (arr: any[]) => any[],
 *   uniform: (a: number, b: number) => number, uint32: () => number, exponential: (rate?: number) => number,
 *   fork: () => object, seed: number }}
 *   next() uniform in [0,1) [-]; normal() standard normal deviate [-]; int(n) integer in 0..n-1 [-];
 *   shuffle(arr) in-place Fisher–Yates shuffle (returns arr); uniform(a,b) uniform in [a,b);
 *   uint32() raw 32-bit output; exponential(rate) exponential deviate; fork() independent child generator.
 */
export function makeRng(seed = 1) {
  if (seed === undefined || seed === null) seed = 1;
  if (typeof seed !== 'number' || !Number.isFinite(seed)) throw new RangeError('makeRng: seed must be a finite number');
  const whole = Math.floor(seed);
  const frac = seed - whole;
  const lo = whole >>> 0;
  const hi = Math.floor(whole / 4294967296) >>> 0;
  const sm = splitmix32((lo ^ Math.imul(hi, 0x85ebca6b) ^ Math.imul((frac * 4294967296) >>> 0, 0xc2b2ae35)) | 0);
  let a = sm() | 0;
  let b = sm() | 0;
  let c = sm() | 0;
  let d = 1;
  const uint32 = () => {
    const t = (((a + b) | 0) + d) | 0;
    d = (d + 1) | 0;
    a = b ^ (b >>> 9);
    b = (c + (c << 3)) | 0;
    c = (c << 21) | (c >>> 11);
    c = (c + t) | 0;
    return t >>> 0;
  };
  for (let i = 0; i < 15; i++) uint32();

  const next = () => {
    const u = uint32() >>> 5; // 27 bits
    const v = uint32() >>> 6; // 26 bits
    return (u * 67108864 + v) / 9007199254740992;
  };

  let spare = 0;
  let hasSpare = false;
  const normal = () => {
    if (hasSpare) {
      hasSpare = false;
      return spare;
    }
    let u;
    let v;
    let s;
    do {
      u = 2 * next() - 1;
      v = 2 * next() - 1;
      s = u * u + v * v;
    } while (s >= 1 || s === 0);
    const m = Math.sqrt((-2 * Math.log(s)) / s);
    spare = v * m;
    hasSpare = true;
    return u * m;
  };

  const int = (n) => {
    if (!(n >= 1)) throw new RangeError('rng.int: n must be >= 1');
    return Math.floor(next() * n);
  };

  const shuffle = (arr) => {
    for (let i = arr.length - 1; i > 0; i--) {
      const j = Math.floor(next() * (i + 1));
      const t = arr[i];
      arr[i] = arr[j];
      arr[j] = t;
    }
    return arr;
  };

  const uniform = (lo2, hi2) => lo2 + (hi2 - lo2) * next();
  const exponential = (rate = 1) => -Math.log(1 - next()) / rate;
  const fork = () => makeRng(uint32() + uint32() / 4294967296);

  return { next, normal, int, shuffle, uniform, uint32, exponential, fork, seed };
}
