// Binomial lattice for options held by the OWNER of an operating asset (expand, abandon).
// The owner receives the asset's cash flow while holding it. With project value V following the risk-neutral
// process dV = (r − q)V dt + σV dW, the cash paid out in a step is V(1 − e^(−qΔt)), and the value of holding is
//   H = V(1 − e^(−qΔt)) + e^(−rΔt) E[F(V', t + Δt)].
// With no options this recursion returns exactly V at every node (V(1 − e^(−qΔt)) + e^(−rΔt) V e^((r−q)Δt) = V),
// so the option value is F(V0, 0) − V0. A lattice that carries only the terminal value, as a claim without
// dividends, understates the value of holding by V(1 − e^(−qT)) and so overstates early exercise when q > 0.
// Ref: Trigeorgis (1996), Real Options, ch. 5; Cox, Ross & Rubinstein (1979) for the lattice parameters.

/**
 * @param {{V0: number, sigma: number, r: number, q: number, T: number, steps: number, american: boolean,
 *   options: ({type: 'expand', factor: number, cost: number} | {type: 'abandon', salvage: number})[]}} o currency,
 *   continuously compounded rates per year, volatility per √year, years. An exercised option ends the others.
 * @returns {{valueWithOptions: number, valueWithoutOptions: number, optionValue: number, exerciseNow: string | null}}
 */
export function assetOptions({ V0, sigma, r, q = 0, T, steps, american = true, options }) {
  if (!(V0 > 0) || !(sigma > 0) || !(T > 0) || !(steps >= 1)) throw new Error('The asset-option lattice needs positive project value, volatility, window and step count.');
  const dt = T / steps; const u = Math.exp(sigma * Math.sqrt(dt)); const d = 1 / u;
  const p = (Math.exp((r - q) * dt) - d) / (u - d);
  if (!(p > 0 && p < 1)) throw new Error('The lattice step is too coarse for this volatility and drift (risk-neutral probability outside 0–1): increase the number of steps.');
  const disc = Math.exp(-r * dt); const pay = 1 - Math.exp(-q * dt);
  const best = (V, hold) => {
    let m = hold; let which = null;
    for (const o of options) {
      const x = o.type === 'expand' ? o.factor * V - o.cost : o.salvage;
      if (x > m) { m = x; which = o.type; }
    }
    return { m, which };
  };
  // node values V0 u^k for k = −steps … steps, tabulated once
  const node = new Float64Array(2 * steps + 1);
  for (let k = -steps; k <= steps; k++) node[k + steps] = V0 * Math.exp(k * sigma * Math.sqrt(dt));
  const val = new Float64Array(steps + 1);
  for (let j = 0; j <= steps; j++) { const V = node[2 * j]; val[j] = best(V, V).m; }
  let exerciseNow = null;
  for (let i = steps - 1; i >= 0; i--) {
    for (let j = 0; j <= i; j++) {
      const V = node[2 * j - i + steps];
      const hold = V * pay + disc * (p * val[j + 1] + (1 - p) * val[j]);
      if (american) { const b = best(V, hold); val[j] = b.m; if (i === 0) exerciseNow = b.which; } else val[j] = hold;
    }
  }
  return { valueWithOptions: val[0], valueWithoutOptions: V0, optionValue: val[0] - V0, exerciseNow };
}
