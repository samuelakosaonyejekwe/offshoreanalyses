/**
 * Generator, converter and drive-train electrical models for offshore renewables (wind, tidal, wave, floating PV).
 *
 * References
 *  - Krause, Wasynczuk, Sudhoff, Analysis of Electric Machinery and Drive Systems (2013): dq models of
 *    permanent-magnet and induction machines (amplitude-invariant transformation, motor sign convention).
 *  - Kloss formula for the torque–slip curve of induction machines (e.g. Boldea & Nasar, The Induction Machine
 *    Handbook, 2002, ch. 7).
 *  - Abad et al., Doubly Fed Induction Machine (2011): steady-state DFIG power relations.
 *  - Infineon AN2008-03 / Semikron Application Manual: averaged conduction and switching losses of a two-level
 *    IGBT bridge with sinusoidal PWM.
 *  - Burton et al., Wind Energy Handbook (2011) §8: optimal-torque (k·ω²) control below rated.
 *  - Esram & Chapman, IEEE Trans. Energy Convers. 22 (2007): perturb-and-observe and incremental-conductance MPPT.
 *
 * Electrical quantities in the dq frame are peak phase values (amplitude-invariant transform), so
 * P = 1.5·(vd·id + vq·iq). Motor convention: positive torque and positive iq drive the shaft; a generator
 * operates with negative iq and negative torque. SI units.
 */

import { cx, cadd, cmul, cdiv, cconj, cabs, cscale, csolve } from './complex.js';

function pos(v, name) {
  if (!(v > 0) || !Number.isFinite(v)) throw new RangeError(`${name}: must be a positive finite number`);
}
function nonNeg(v, name) {
  if (!(v >= 0) || !Number.isFinite(v)) throw new RangeError(`${name}: must be a non-negative finite number`);
}
function fin(v, name) {
  if (!Number.isFinite(v)) throw new RangeError(`${name}: must be finite`);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Permanent-magnet synchronous generator                                                                       */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Electromagnetic torque of a PMSM: Te = 1.5·p·(ψm·iq + (Ld − Lq)·id·iq).
 * @param {object} o
 * @param {number} o.polePairs p @param {number} o.psiM permanent-magnet flux linkage (peak), Wb
 * @param {number} o.Ld d-axis inductance, H @param {number} o.Lq q-axis inductance, H
 * @param {number} o.id d-axis current (peak), A @param {number} o.iq q-axis current (peak), A
 * @returns {number} torque, N·m
 */
export function pmsmTorque({ polePairs, psiM, Ld, Lq, id, iq }) {
  pos(polePairs, 'polePairs'); nonNeg(psiM, 'psiM'); pos(Ld, 'Ld'); pos(Lq, 'Lq'); fin(id, 'id'); fin(iq, 'iq');
  return 1.5 * polePairs * (psiM * iq + (Ld - Lq) * id * iq);
}

/**
 * Current references for a torque demand.
 *  - 'zeroD': id = 0, iq = Te/(1.5·p·ψm).
 *  - 'mtpa' (maximum torque per ampere): for a current magnitude Is,
 *    id = [−ψm + sqrt(ψm² + 8·(Ld − Lq)²·Is²)]/(4·(Ld − Lq)), iq = ±sqrt(Is² − id²); Is is found by bisection so
 *    that the torque equals the demand. For Ld = Lq it reduces to zero-d-axis control.
 * Reference: Krause et al.; Morimoto et al., IEEE Trans. Ind. Appl. 26 (1990) for the MTPA expression.
 * @param {object} o machine parameters polePairs, psiM, Ld, Lq and
 * @param {number} o.torque demanded electromagnetic torque (negative when generating), N·m
 * @param {'zeroD'|'mtpa'} [o.strategy='mtpa']
 * @returns {{id:number, iq:number, Is:number}} peak currents, A
 */
export function pmsmCurrentReference({ polePairs, psiM, Ld, Lq, torque, strategy = 'mtpa' }) {
  pos(polePairs, 'polePairs'); pos(psiM, 'psiM'); pos(Ld, 'Ld'); pos(Lq, 'Lq'); fin(torque, 'torque');
  const sign = torque < 0 ? -1 : 1;
  const T = Math.abs(torque);
  if (strategy === 'zeroD' || Math.abs(Ld - Lq) < 1e-12 * Math.max(Ld, Lq)) {
    if (strategy !== 'zeroD' && strategy !== 'mtpa') throw new RangeError('strategy: must be zeroD or mtpa');
    const iq = T / (1.5 * polePairs * psiM);
    return { id: 0, iq: sign * iq, Is: iq };
  }
  if (strategy !== 'mtpa') throw new RangeError('strategy: must be zeroD or mtpa');
  const dL = Ld - Lq;
  const at = (Is) => {
    const id = (-psiM + Math.sqrt(psiM * psiM + 8 * dL * dL * Is * Is)) / (4 * dL);
    const iq = Math.sqrt(Math.max(0, Is * Is - id * id));
    return { id, iq, T: 1.5 * polePairs * (psiM * iq + dL * id * iq) };
  };
  let lo = 0; let hi = T / (1.5 * polePairs * psiM) + 1e-9;
  while (at(hi).T < T) hi *= 2;
  for (let i = 0; i < 100; i++) { const mid = 0.5 * (lo + hi); if (at(mid).T < T) lo = mid; else hi = mid; }
  const r = at(0.5 * (lo + hi));
  return { id: r.id, iq: sign * r.iq, Is: 0.5 * (lo + hi) };
}

/**
 * Steady-state operating point of a permanent-magnet synchronous machine in the rotor dq frame:
 *   vd = Rs·id − ωe·Lq·iq,  vq = Rs·iq + ωe·(Ld·id + ψm),  ωe = p·ωm,
 *   Te = 1.5·p·(ψm·iq + (Ld − Lq)·id·iq),  Pe = 1.5·(vd·id + vq·iq)  (terminal power, negative = generated).
 * Losses: copper 1.5·Rs·(id² + iq²); iron kh·ωe + ke·ωe² scaled by (ψ/ψm)² with ψ the stator flux magnitude
 * (hysteresis and eddy-current coefficients in W·s/rad and W·s²/rad² at rated flux, from the machine data sheet
 * or a no-load test); mechanical friction/windage Tf0·ωm + kw·ωm³. The shaft power is
 * Pshaft = Te·ωm − iron − mechanical (motor convention), so for a generator |Pshaft| = |Te·ωm| + iron + mechanical.
 * Reference: Krause et al. ch. 4 (motor convention, amplitude-invariant dq).
 * @param {object} o
 * @param {number} o.polePairs @param {number} o.psiM Wb @param {number} o.Ld H @param {number} o.Lq H
 * @param {number} o.Rs stator phase resistance, Ω
 * @param {number} o.speed mechanical speed ωm, rad/s
 * @param {number} [o.torque] electromagnetic torque demand, N·m (negative = generating); or give id and iq
 * @param {number} [o.id] @param {number} [o.iq] peak dq currents, A
 * @param {'zeroD'|'mtpa'} [o.strategy='mtpa']
 * @param {number} [o.kh=0] @param {number} [o.ke=0] iron-loss coefficients at rated flux
 * @param {number} [o.frictionTorque=0] N·m @param {number} [o.windage=0] W·s³/rad³
 * @returns {{id:number, iq:number, vd:number, vq:number, Vs:number, Is:number, Vll:number, Irms:number,
 *            torque:number, electricalPower:number, shaftPower:number, copperLoss:number, ironLoss:number,
 *            mechanicalLoss:number, efficiency:number, powerFactor:number, frequency:number}}
 *          peak dq values (A, V), peak phase voltage/current, rms line-to-line voltage (V) and rms phase current
 *          (A), torque (N·m), powers and losses (W), efficiency (output/input in the direction of power flow),
 *          terminal power factor and electrical frequency (Hz)
 */
export function pmsmSteadyState({ polePairs, psiM, Ld, Lq, Rs, speed, torque, id, iq, strategy = 'mtpa', kh = 0, ke = 0,
  frictionTorque = 0, windage = 0 }) {
  pos(polePairs, 'polePairs'); pos(psiM, 'psiM'); pos(Ld, 'Ld'); pos(Lq, 'Lq'); nonNeg(Rs, 'Rs'); nonNeg(speed, 'speed');
  nonNeg(kh, 'kh'); nonNeg(ke, 'ke'); nonNeg(frictionTorque, 'frictionTorque'); nonNeg(windage, 'windage');
  let d = id; let q = iq;
  if (d === undefined || q === undefined) {
    if (torque === undefined) throw new RangeError('torque: required when id and iq are not given');
    ({ id: d, iq: q } = pmsmCurrentReference({ polePairs, psiM, Ld, Lq, torque, strategy }));
  }
  const we = polePairs * speed;
  const vd = Rs * d - we * Lq * q;
  const vq = Rs * q + we * (Ld * d + psiM);
  const Te = 1.5 * polePairs * (psiM * q + (Ld - Lq) * d * q);
  const Pe = 1.5 * (vd * d + vq * q);
  const copperLoss = 1.5 * Rs * (d * d + q * q);
  const flux = Math.hypot(Ld * d + psiM, Lq * q);
  const ironLoss = (kh * we + ke * we * we) * (flux / psiM) ** 2;
  const mechanicalLoss = frictionTorque * speed + windage * speed ** 3;
  const shaftPower = Te * speed - ironLoss - mechanicalLoss;
  const Vs = Math.hypot(vd, vq); const Is = Math.hypot(d, q);
  // Generating: shaft input |shaftPower|, electrical output −Pe. Motoring: electrical input Pe (+ iron drawn
  // from the supply is already inside the model as a shaft-side loss), shaft output shaftPower.
  let efficiency = NaN;
  if (Pe < 0 && shaftPower < 0) efficiency = -Pe / -shaftPower;
  else if (Pe > 0 && shaftPower > 0) efficiency = shaftPower / Pe;
  return { id: d, iq: q, vd, vq, Vs, Is, Vll: Vs * Math.sqrt(1.5), Irms: Is / Math.SQRT2, torque: Te, electricalPower: Pe,
    shaftPower, copperLoss, ironLoss, mechanicalLoss, efficiency,
    powerFactor: Vs * Is > 0 ? Pe / (1.5 * Vs * Is) : NaN, frequency: we / (2 * Math.PI) };
}

/**
 * Generator efficiency map over speed and torque (generating quadrant).
 * @param {object} o machine and loss parameters of pmsmSteadyState plus
 * @param {number[]} o.speeds mechanical speeds, rad/s @param {number[]} o.torques shaft-torque magnitudes, N·m
 * @returns {{speeds:number[], torques:number[], efficiency:number[][], power:number[][]}} efficiency[iTorque][iSpeed]
 *          (NaN where no power is generated) and electrical output power (W)
 */
export function pmsmEfficiencyMap({ speeds, torques, ...machine }) {
  if (!Array.isArray(speeds) || !Array.isArray(torques)) throw new RangeError('speeds, torques: must be arrays');
  const efficiency = []; const power = [];
  for (const T of torques) {
    const rowE = []; const rowP = [];
    for (const w of speeds) {
      const r = pmsmSteadyState({ ...machine, speed: w, torque: -Math.abs(T) });
      rowE.push(r.efficiency); rowP.push(-r.electricalPower);
    }
    efficiency.push(rowE); power.push(rowP);
  }
  return { speeds: [...speeds], torques: [...torques], efficiency, power };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Induction machines                                                                                           */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Torque–slip characteristic of a squirrel-cage induction machine from the equivalent circuit (stator referred):
 *   T(s) = 3·p/ωs · Vth²·(R2/s) / [(Rth + R2/s)² + (Xth + X2)²],
 * with the Thévenin source of the stator and magnetising branch. Negative slip gives generator operation.
 * The Kloss form T/Tb = (2 + 2·a·sb)/(s/sb + sb/s + 2·a·sb), a = Rth/R2, Tb the motoring breakdown torque at
 * sb = R2/sqrt(Rth² + (Xth + X2)²), is algebraically identical to the circuit expression and is returned as a
 * cross-check; a = 0 gives the simplified form 2/(s/sb + sb/s) of klossTorque.
 * Reference: Boldea & Nasar ch. 7; Fitzgerald, Kingsley & Umans, Electric Machinery, ch. 6.
 * @param {object} o
 * @param {number} o.V line-to-line supply voltage (rms), V @param {number} o.f supply frequency, Hz
 * @param {number} o.polePairs @param {number} o.R1 stator resistance, Ω @param {number} o.X1 stator leakage, Ω
 * @param {number} o.R2 rotor resistance (referred), Ω @param {number} o.X2 rotor leakage (referred), Ω
 * @param {number} o.Xm magnetising reactance, Ω
 * @param {number} o.slip slip s = (ωs − p·ωm)/ωs (non-zero; negative = generating)
 * @returns {{torque:number, torqueKloss:number, breakdownTorque:number, breakdownSlip:number, I1:number,
 *            powerFactor:number, P:number, Q:number, speed:number, efficiency:number}} N·m, A (rms), stator
 *          electrical power (W, positive = absorbed) and reactive power (var, positive = absorbed), mechanical
 *          speed (rad/s), efficiency in the direction of power flow (iron and friction losses not included)
 */
export function inductionMachine({ V, f, polePairs, R1, X1, R2, X2, Xm, slip }) {
  pos(V, 'V'); pos(f, 'f'); pos(polePairs, 'polePairs'); nonNeg(R1, 'R1'); nonNeg(X1, 'X1'); pos(R2, 'R2'); nonNeg(X2, 'X2'); pos(Xm, 'Xm');
  if (!Number.isFinite(slip) || slip === 0) throw new RangeError('slip: must be finite and non-zero');
  const ws = (2 * Math.PI * f) / polePairs; // synchronous mechanical speed
  const V1 = V / Math.sqrt(3);
  const z1 = cx(R1, X1); const zm = cx(0, Xm);
  const vth = cabs(cdiv(cmul(cx(V1, 0), zm), cadd(z1, zm)));
  const zth = cdiv(cmul(z1, zm), cadd(z1, zm));
  const xs = zth.im + X2;
  const tq = (s) => ((3 / ws) * vth * vth * (R2 / s)) / ((zth.re + R2 / s) ** 2 + xs * xs);
  const sb = R2 / Math.hypot(zth.re, xs);
  // Breakdown torque for the sign of the operating slip (motoring and generating maxima differ when Rth > 0).
  const sgn = slip > 0 ? 1 : -1;
  const Tb = tq(sgn * sb);
  const a = zth.re / R2;
  const torqueKloss = (tq(sb) * (2 + 2 * a * sb)) / (slip / sb + sb / slip + 2 * a * sb);
  const z2 = cx(R2 / slip, X2);
  const zin = cadd(z1, cdiv(cmul(zm, z2), cadd(zm, z2)));
  const I1 = cdiv(cx(V1, 0), zin);
  const S = cscale(cmul(cx(V1, 0), cconj(I1)), 3);
  const torque = tq(slip);
  const speed = ws * (1 - slip);
  const Pmech = torque * speed;
  let efficiency = NaN;
  if (S.re > 0 && Pmech > 0) efficiency = Pmech / S.re; else if (S.re < 0 && Pmech < 0) efficiency = S.re / Pmech;
  return { torque, torqueKloss, breakdownTorque: Tb, breakdownSlip: sgn * sb, I1: cabs(I1), powerFactor: S.re / cabs(S),
    P: S.re, Q: S.im, speed, efficiency };
}

/**
 * Kloss torque–slip relation (stator resistance neglected): T/Tb = 2/(s/sb + sb/s).
 * @param {object} o
 * @param {number} o.slip @param {number} o.breakdownSlip @param {number} o.breakdownTorque N·m
 * @returns {number} torque, N·m (sign follows the slip)
 */
export function klossTorque({ slip, breakdownSlip, breakdownTorque }) {
  if (!Number.isFinite(slip) || slip === 0) throw new RangeError('slip: must be finite and non-zero');
  pos(breakdownSlip, 'breakdownSlip'); pos(breakdownTorque, 'breakdownTorque');
  return (2 * breakdownTorque) / (slip / breakdownSlip + breakdownSlip / slip);
}

/**
 * Doubly-fed induction generator in steady state. Mesh equations of the T equivalent circuit (motor convention,
 * rotor referred to the stator, rotor circuit divided by slip):
 *   Vs = (Rs + jXs)·Is + jXm·(Is + Ir),   Vr/s = (Rr/s + jXr)·Ir + jXm·(Is + Ir).
 * Powers: Ps = 3·Re(Vs·Is*), Pr = 3·Re(Vr·Ir*) (both positive when absorbed) and Pmech = Ps + Pr − copper losses
 * (motor convention). Neglecting losses Pr = −s·Ps, so the net electrical power is (1 − s)·Ps: a generator
 * (Ps < 0) absorbs rotor power below synchronous speed (s > 0) and delivers it above (s < 0).
 * Reference: Abad et al. (2011) ch. 3.
 * @param {object} o
 * @param {number} o.Vs stator line-to-line voltage (rms), V
 * @param {{re:number, im:number}} o.Vr rotor phase voltage phasor referred to the stator (rms, stator voltage
 *        as phase reference), V
 * @param {number} o.slip (non-zero) @param {number} o.Rs Ω @param {number} o.Xs stator leakage reactance, Ω
 * @param {number} o.Rr Ω @param {number} o.Xr rotor leakage reactance at stator frequency, Ω @param {number} o.Xm Ω
 * @returns {{Is:{re:number, im:number}, Ir:{re:number, im:number}, Ps:number, Qs:number, Pr:number, Qr:number,
 *            copperLoss:number, Pmech:number, Ptotal:number, rotorPowerRatio:number}} A (rms phasors), W, var;
 *          Pmech is the mechanical power delivered to the shaft (negative when driven as a generator), Ptotal =
 *          Ps + Pr the net electrical power absorbed, rotorPowerRatio = Pr/Ps (≈ −s when losses are small)
 */
export function dfigSteadyState({ Vs, Vr, slip, Rs, Xs, Rr, Xr, Xm }) {
  pos(Vs, 'Vs'); nonNeg(Rs, 'Rs'); nonNeg(Xs, 'Xs'); nonNeg(Rr, 'Rr'); nonNeg(Xr, 'Xr'); pos(Xm, 'Xm');
  if (!Number.isFinite(slip) || slip === 0) throw new RangeError('slip: must be finite and non-zero');
  const V1 = cx(Vs / Math.sqrt(3), 0);
  const A = [[cx(Rs, Xs + Xm), cx(0, Xm)], [cx(0, Xm), cx(Rr / slip, Xr + Xm)]];
  const [Is, Ir] = csolve(A, [V1, cscale(Vr, 1 / slip)]);
  const Ss = cscale(cmul(V1, cconj(Is)), 3);
  const Sr = cscale(cmul(Vr, cconj(Ir)), 3);
  const copperLoss = 3 * (Rs * cabs(Is) ** 2 + Rr * cabs(Ir) ** 2);
  const Ptotal = Ss.re + Sr.re;
  return { Is, Ir, Ps: Ss.re, Qs: Ss.im, Pr: Sr.re, Qr: Sr.im, copperLoss, Pmech: Ptotal - copperLoss, Ptotal,
    rotorPowerRatio: Ss.re !== 0 ? Sr.re / Ss.re : NaN };
}

/**
 * Wound-field (or any round-rotor/salient) synchronous machine, steady state, resistance neglected:
 * P = 3·E·V·sin δ/Xd + 3·V²/2·(1/Xq − 1/Xd)·sin 2δ,  Q = 3·E·V·cos δ/Xd − 3·V²·(cos²δ/Xd + sin²δ/Xq)
 * (generator convention: P, Q delivered to the terminals; E, V phase rms). For Xq = Xd the round-rotor
 * expressions result.
 * Reference: Fitzgerald, Kingsley & Umans ch. 5.
 * @param {object} o
 * @param {number} o.E internal (excitation) voltage, line-to-line rms, V @param {number} o.V terminal voltage,
 *        line-to-line rms, V @param {number} o.delta load angle, rad @param {number} o.Xd Ω @param {number} [o.Xq=Xd] Ω
 * @returns {{P:number, Q:number, Pmax:number}} W, var; Pmax is the round-rotor pull-out power 3·E·V/Xd (phase values)
 */
export function synchronousMachinePower({ E, V, delta, Xd, Xq }) {
  nonNeg(E, 'E'); pos(V, 'V'); fin(delta, 'delta'); pos(Xd, 'Xd');
  const xq = Xq ?? Xd; pos(xq, 'Xq');
  const e = E / Math.sqrt(3); const v = V / Math.sqrt(3);
  const P = (3 * e * v * Math.sin(delta)) / Xd + 1.5 * v * v * (1 / xq - 1 / Xd) * Math.sin(2 * delta);
  const Q = (3 * e * v * Math.cos(delta)) / Xd - 3 * v * v * (Math.cos(delta) ** 2 / Xd + Math.sin(delta) ** 2 / xq);
  return { P, Q, Pmax: (3 * e * v) / Xd };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Converters                                                                                                   */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Averaged losses of a three-phase two-level IGBT bridge with sinusoidal PWM (six switches and six diodes):
 *   IGBT conduction (each)  = Vce0·Î·(1/(2π) + m·cosφ/8) + rce·Î²·(1/8 + m·cosφ/(3π))
 *   diode conduction (each) = Vf0·Î·(1/(2π) − m·cosφ/8) + rf·Î²·(1/8 − m·cosφ/(3π))
 *   switching (each pair)   = fsw·(Eon + Eoff + Erec)/π · (Î/Iref)·(Vdc/Vref)
 * Î is the peak phase current, m the modulation index, cosφ the displacement factor seen by the bridge
 * (negative when the bridge rectifies). Device parameters come from the module data sheet; the defaults describe
 * a generic 1700 V / 1000 A-class module and are placeholders to be edited, not data for any product.
 * The result is also returned as the polynomial P = a1·Î + a2·Î² in the peak current.
 * Reference: Infineon AN2008-03; Semikron Application Manual Power Semiconductors §5.2.
 * @param {object} o
 * @param {number} o.Irms phase current (rms), A @param {number} o.Vdc DC-link voltage, V
 * @param {number} [o.modulationIndex=0.9] @param {number} [o.cosPhi=1] @param {number} [o.fsw=2500] Hz
 * @param {number} [o.Vce0=1.0] V @param {number} [o.rce=1.2e-3] Ω @param {number} [o.Vf0=1.0] V @param {number} [o.rf=0.8e-3] Ω
 * @param {number} [o.Eon=0.35] J @param {number} [o.Eoff=0.40] J @param {number} [o.Erec=0.25] J at Iref, Vref
 * @param {number} [o.Iref=1000] A @param {number} [o.Vref=900] V @param {number} [o.parallel=1] modules in parallel
 * @returns {{conduction:number, switching:number, total:number, a1:number, a2:number}} W; a1 in W/A, a2 in W/A²
 */
export function converterLoss({ Irms, Vdc, modulationIndex = 0.9, cosPhi = 1, fsw = 2500, Vce0 = 1.0, rce = 1.2e-3, Vf0 = 1.0,
  rf = 0.8e-3, Eon = 0.35, Eoff = 0.40, Erec = 0.25, Iref = 1000, Vref = 900, parallel = 1 }) {
  nonNeg(Irms, 'Irms'); pos(Vdc, 'Vdc'); nonNeg(fsw, 'fsw'); pos(Iref, 'Iref'); pos(Vref, 'Vref'); pos(parallel, 'parallel');
  if (!(modulationIndex >= 0 && modulationIndex <= 1.1547)) throw new RangeError('modulationIndex: must be in [0, 1.155]');
  if (!(cosPhi >= -1 && cosPhi <= 1)) throw new RangeError('cosPhi: must be in [−1, 1]');
  const mc = modulationIndex * cosPhi;
  // Per-device coefficients; `parallel` modules share the current equally.
  const lin = Vce0 * (1 / (2 * Math.PI) + mc / 8) + Vf0 * (1 / (2 * Math.PI) - mc / 8)
    + ((fsw * (Eon + Eoff + Erec)) / Math.PI) * (1 / Iref) * (Vdc / Vref);
  const quad = (rce * (1 / 8 + mc / (3 * Math.PI)) + rf * (1 / 8 - mc / (3 * Math.PI))) / parallel;
  const a1 = 6 * lin; const a2 = 6 * quad;
  const Ipk = Irms * Math.SQRT2;
  const switching = 6 * ((fsw * (Eon + Eoff + Erec)) / Math.PI) * (Ipk / Iref) * (Vdc / Vref);
  const total = a1 * Ipk + a2 * Ipk * Ipk;
  return { conduction: total - switching, switching, total, a1, a2 };
}

/**
 * DC-link capacitor sizing rules.
 *  - Energy (inertia) rule: C = 2·S·τ/Vdc², τ = stored energy / rated power (typically 5–30 ms for wind-turbine
 *    back-to-back converters).
 *  - Ripple rule for a power imbalance ΔP at angular frequency ω (e.g. twice the grid frequency for an
 *    unbalanced fault, or the machine-side pulsation): C = ΔP/(ω·Vdc·ΔVdc), ΔVdc peak-to-peak.
 * The larger governs. Reference: Mohan, Undeland & Robbins, Power Electronics (2003) §8.
 * @param {object} o
 * @param {number} o.S rated power, VA @param {number} o.Vdc V @param {number} [o.tau=0.01] s
 * @param {number} [o.rippleRatio] ΔP as a fraction of S @param {number} [o.rippleFrequency] Hz
 * @param {number} [o.rippleVoltageFraction=0.05] ΔVdc/Vdc (peak-to-peak)
 * @returns {{Cenergy:number, Cripple:number, C:number, storedEnergy:number}} F, F, F, J
 */
export function dcLinkCapacitance({ S, Vdc, tau = 0.01, rippleRatio, rippleFrequency, rippleVoltageFraction = 0.05 }) {
  pos(S, 'S'); pos(Vdc, 'Vdc'); pos(tau, 'tau'); pos(rippleVoltageFraction, 'rippleVoltageFraction');
  const Cenergy = (2 * S * tau) / (Vdc * Vdc);
  let Cripple = 0;
  if (rippleRatio !== undefined && rippleFrequency !== undefined) {
    nonNeg(rippleRatio, 'rippleRatio'); pos(rippleFrequency, 'rippleFrequency');
    Cripple = (rippleRatio * S) / (2 * Math.PI * rippleFrequency * Vdc * Vdc * rippleVoltageFraction);
  }
  const C = Math.max(Cenergy, Cripple);
  return { Cenergy, Cripple, C, storedEnergy: 0.5 * C * Vdc * Vdc };
}

/**
 * Reactive-power capability (PQ diagram) of a grid-side converter behind a reactance X:
 *  - current limit: P² + Q² ≤ (√3·V·Imax)²;
 *  - converter voltage limit: P² + (Q + V²/X)² ≤ (V·Vc,max/X)² (circle centred at (0, −V²/X)), Vc,max being the
 *    largest fundamental line-to-line voltage the converter can synthesise (m,max·Vdc·√3/(2√2) for sinusoidal PWM).
 * Q > 0 is reactive power exported (over-excited). Returns the limits at a given P.
 * Reference: e.g. Engelhardt et al., IEEE Trans. Energy Convers. 26 (2011) for wind-turbine capability curves.
 * @param {object} o
 * @param {number} o.V grid line-to-line voltage, V @param {number} o.Imax converter current limit (rms), A
 * @param {number} o.P active power, W @param {number} [o.X] coupling reactance (filter + transformer), Ω
 * @param {number} [o.VcMax] maximum converter line-to-line voltage (rms), V
 * @returns {{Smax:number, Qmax:number, Qmin:number, limitedBy:string, feasible:boolean}} VA, var
 */
export function converterCapability({ V, Imax, P, X, VcMax }) {
  pos(V, 'V'); pos(Imax, 'Imax'); fin(P, 'P');
  const Smax = Math.sqrt(3) * V * Imax;
  if (Math.abs(P) > Smax) return { Smax, Qmax: 0, Qmin: 0, limitedBy: 'current', feasible: false };
  const qI = Math.sqrt(Smax * Smax - P * P);
  let Qmax = qI; let limitedBy = 'current';
  if (X !== undefined && VcMax !== undefined) {
    pos(X, 'X'); pos(VcMax, 'VcMax');
    const rad = (V * VcMax) / X;
    if (Math.abs(P) > rad) return { Smax, Qmax: 0, Qmin: -qI, limitedBy: 'voltage', feasible: false };
    const qV = Math.sqrt(rad * rad - P * P) - (V * V) / X;
    if (qV < Qmax) { Qmax = qV; limitedBy = 'voltage'; }
  }
  return { Smax, Qmax, Qmin: -qI, limitedBy, feasible: Qmax >= -qI };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Drive train and operating schedule                                                                           */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Gearbox efficiency from stage losses: each stage loses a fraction of the transmitted power (about 1 % per
 * planetary stage and 1.5–2 % per parallel-shaft stage at rated load — arguments), plus a constant no-load loss
 * as a fraction of rated power (churning, seals, bearings). Direct drive: no stages, efficiency 1.
 * η(P) = [P·Π(1 − l_i) − p0·Prated]/P.
 * Reference: loss breakdown practice, e.g. IEC 61400-4 informative guidance; coefficients are user data.
 * @param {object} o
 * @param {number[]} [o.stageLosses=[]] load-dependent loss fraction per stage
 * @param {number} [o.noLoadLossFraction=0] constant loss / rated power
 * @param {number} o.power input shaft power, W @param {number} o.ratedPower W
 * @returns {{efficiency:number, outputPower:number, loss:number}}
 */
export function gearboxEfficiency({ stageLosses = [], noLoadLossFraction = 0, power, ratedPower }) {
  nonNeg(power, 'power'); pos(ratedPower, 'ratedPower'); nonNeg(noLoadLossFraction, 'noLoadLossFraction');
  let eta = 1;
  for (const l of stageLosses) { if (!(l >= 0 && l < 1)) throw new RangeError('stageLosses: fractions must be in [0, 1)'); eta *= 1 - l; }
  const outputPower = Math.max(0, power * eta - noLoadLossFraction * ratedPower);
  return { efficiency: power > 0 ? outputPower / power : 0, outputPower, loss: power - outputPower };
}

/**
 * Optimal-torque (maximum power point tracking) law for a wind or tidal turbine below rated:
 * T = k·ω², k = ½·ρ·π·R⁵·Cp,max/λopt³ (rotor side); on the generator side of a gearbox of ratio N, k/N³.
 * Reference: Burton et al. §8.3.
 * @param {object} o
 * @param {number} o.rho fluid density, kg/m³ @param {number} o.R rotor radius, m
 * @param {number} o.CpMax maximum power coefficient @param {number} o.lambdaOpt optimal tip-speed ratio
 * @param {number} [o.gearRatio=1] generator speed / rotor speed
 * @param {number} [o.omega] rotor speed for which the torque and power are wanted, rad/s
 * @returns {{k:number, kGenerator:number, torque:?number, power:?number}} k in N·m·s²/rad²; rotor torque (N·m) and
 *          power (W) at omega
 */
export function optimalTorqueGain({ rho, R, CpMax, lambdaOpt, gearRatio = 1, omega }) {
  pos(rho, 'rho'); pos(R, 'R'); pos(CpMax, 'CpMax'); pos(lambdaOpt, 'lambdaOpt'); pos(gearRatio, 'gearRatio');
  const k = (0.5 * rho * Math.PI * R ** 5 * CpMax) / lambdaOpt ** 3;
  return { k, kGenerator: k / gearRatio ** 3, torque: omega !== undefined ? k * omega * omega : null,
    power: omega !== undefined ? k * omega ** 3 : null };
}

/**
 * Electrical power and efficiency chain over a variable-speed operating schedule: for each pair (rotor speed,
 * aerodynamic/hydrodynamic torque) the shaft power passes through the gearbox, the PMSG (dq model with losses),
 * the back-to-back converter (two bridges, averaged loss model) and the transformer.
 * @param {object} o
 * @param {number[]} o.rotorSpeed rad/s @param {number[]} o.rotorTorque N·m (magnitudes)
 * @param {object} o.generator parameters of pmsmSteadyState (polePairs, psiM, Ld, Lq, Rs, kh, ke, ...)
 * @param {{ratio?:number, stageLosses?:number[], noLoadLossFraction?:number}} [o.gearbox] ratio 1 = direct drive
 * @param {object} [o.converter] parameters of converterLoss (Vdc required) — omitted: lossless converter
 * @param {number} [o.gridVoltage] line-to-line voltage at the grid-side bridge, V (for its current; default: same
 *        current as the generator side)
 * @param {{Sr:number, Pk:number, P0:number}} [o.transformer] rated power (VA), load and no-load loss (W)
 * @param {number} o.ratedPower rated shaft power, W
 * @returns {{shaftPower:number[], generatorPower:number[], gridPower:number[], efficiency:number[],
 *            losses:{gearbox:number[], generator:number[], converter:number[], transformer:number[]}}} W
 */
export function operatingSchedule({ rotorSpeed, rotorTorque, generator, gearbox = {}, converter, gridVoltage, transformer, ratedPower }) {
  if (!Array.isArray(rotorSpeed) || !Array.isArray(rotorTorque) || rotorSpeed.length !== rotorTorque.length) {
    throw new RangeError('rotorSpeed, rotorTorque: must be arrays of equal length');
  }
  pos(ratedPower, 'ratedPower');
  const ratio = gearbox.ratio ?? 1;
  const out = { shaftPower: [], generatorPower: [], gridPower: [], efficiency: [],
    losses: { gearbox: [], generator: [], converter: [], transformer: [] } };
  for (let i = 0; i < rotorSpeed.length; i++) {
    const w = rotorSpeed[i]; const T = Math.abs(rotorTorque[i]);
    const Pshaft = w * T;
    const gb = gearboxEfficiency({ stageLosses: gearbox.stageLosses ?? [], noLoadLossFraction: gearbox.noLoadLossFraction ?? 0,
      power: Pshaft, ratedPower });
    const wg = w * ratio;
    let Pgen = 0; let genLoss = 0; let Irms = 0;
    if (gb.outputPower > 0 && wg > 0) {
      // Find the electromagnetic torque whose shaft-side power equals the gearbox output (iron and friction
      // losses are on the shaft side of the air gap).
      const target = gb.outputPower;
      const shaftIn = (Te) => -pmsmSteadyState({ ...generator, speed: wg, torque: -Te }).shaftPower;
      let lo = 0; let hi = target / wg;
      if (shaftIn(0) >= target) { lo = 0; hi = 0; } else {
        for (let k = 0; k < 60; k++) { const mid = 0.5 * (lo + hi); if (shaftIn(mid) < target) lo = mid; else hi = mid; }
      }
      const st = pmsmSteadyState({ ...generator, speed: wg, torque: -0.5 * (lo + hi) });
      Pgen = Math.max(0, -st.electricalPower);
      genLoss = gb.outputPower - Pgen;
      Irms = st.Irms;
    } else genLoss = gb.outputPower;
    let convLoss = 0;
    if (converter && Pgen > 0) {
      const machineSide = converterLoss({ ...converter, Irms, cosPhi: -Math.abs(converter.cosPhiMachine ?? 1) }).total;
      const Igrid = gridVoltage ? Pgen / (Math.sqrt(3) * gridVoltage) : Irms;
      const gridSide = converterLoss({ ...converter, Irms: Igrid, cosPhi: converter.cosPhiGrid ?? 1 }).total;
      convLoss = Math.min(Pgen, machineSide + gridSide);
    }
    const Pconv = Pgen - convLoss;
    let trLoss = 0;
    if (transformer && Pconv > 0) {
      const x = Pconv / transformer.Sr;
      trLoss = Math.min(Pconv, transformer.P0 + x * x * transformer.Pk);
    }
    const Pgrid = Pconv - trLoss;
    out.shaftPower.push(Pshaft); out.generatorPower.push(Pgen); out.gridPower.push(Pgrid);
    out.efficiency.push(Pshaft > 0 ? Pgrid / Pshaft : 0);
    out.losses.gearbox.push(gb.loss); out.losses.generator.push(genLoss); out.losses.converter.push(convLoss);
    out.losses.transformer.push(trLoss);
  }
  return out;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Maximum power point tracking for PV                                                                          */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Discrete simulation of a PV maximum-power-point tracker on a current–voltage characteristic supplied by the
 * caller.
 *  - 'perturbObserve': the voltage reference is stepped by ±ΔV; the direction is reversed when the power fell.
 *  - 'incrementalConductance': the step direction follows the sign of dI/dV + I/V (zero at the maximum power
 *    point), and no step is taken when it is within the tolerance.
 * Reference: Esram & Chapman (2007).
 * @param {object} o
 * @param {(V:number, step:number)=>number} o.current array current as a function of voltage (A); `step` lets the
 *        caller change irradiance during the run
 * @param {number} o.V0 initial voltage, V @param {number} o.dV perturbation step, V @param {number} [o.steps=200]
 * @param {'perturbObserve'|'incrementalConductance'} [o.method='perturbObserve']
 * @param {number} [o.tolerance=1e-3] relative tolerance on |dI/dV + I/V|·V/I for incremental conductance
 * @param {number} [o.Vmin=0] @param {number} [o.Vmax=Infinity] voltage window of the converter, V
 * @returns {{V:Float64Array, P:Float64Array, Vfinal:number, Pfinal:number, Pmean:number}} trajectories and the
 *          mean power over the last quarter of the run (W), which includes the steady-state oscillation loss
 */
export function pvMppt({ current, V0, dV, steps = 200, method = 'perturbObserve', tolerance = 1e-3, Vmin = 0, Vmax = Infinity }) {
  if (typeof current !== 'function') throw new RangeError('current: must be a function I(V, step)');
  pos(V0, 'V0'); pos(dV, 'dV');
  if (!(steps >= 2)) throw new RangeError('steps: at least 2');
  if (method !== 'perturbObserve' && method !== 'incrementalConductance') throw new RangeError('method: must be perturbObserve or incrementalConductance');
  const V = new Float64Array(steps); const P = new Float64Array(steps);
  let v = V0; let dir = 1; let vPrev = V0; let iPrev = current(V0, 0); let pPrev = vPrev * iPrev;
  V[0] = v; P[0] = pPrev;
  v = Math.min(Vmax, Math.max(Vmin, v + dV));
  for (let k = 1; k < steps; k++) {
    const i = current(v, k);
    const p = v * i;
    V[k] = v; P[k] = p;
    let move;
    if (method === 'perturbObserve') {
      if (p < pPrev) dir = -dir;
      move = dir * dV;
    } else {
      const dv = v - vPrev; const di = i - iPrev;
      if (dv === 0) move = di === 0 ? 0 : Math.sign(di) * dV;
      else {
        const g = di / dv + i / v; // dP/dV divided by V
        move = Math.abs(g) * v <= tolerance * Math.abs(i) ? 0 : Math.sign(g) * dV;
      }
    }
    vPrev = v; iPrev = i; pPrev = p;
    v = Math.min(Vmax, Math.max(Vmin, v + move));
  }
  let sum = 0; const from = Math.floor(0.75 * steps);
  for (let k = from; k < steps; k++) sum += P[k];
  return { V, P, Vfinal: V[steps - 1], Pfinal: P[steps - 1], Pmean: sum / (steps - from) };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Wave-energy power take-off                                                                                   */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Permanent-magnet linear generator feeding a resistive load (per-phase equivalent, synchronous reactance
 * neglected for the slow motion of wave devices): EMF e = ke·v (peak phase), force constant kf = 1.5·ke for a
 * three-phase machine in the amplitude-invariant frame, current i = e/(Rg + RL), electromagnetic damping
 * c = F/v = 1.5·ke²/(Rg + RL), electrical output 1.5·i²·RL and efficiency RL/(Rg + RL).
 * For sinusoidal motion of velocity amplitude v̂ the mean output is half the peak value.
 * Reference: Polinder et al., IEEE Trans. Energy Convers. 19 (2004) (linear PM generator for the AWS), reduced to
 * its resistive-load limit.
 * @param {object} o
 * @param {number} o.ke EMF constant, V·s/m (peak phase voltage per unit translator speed)
 * @param {number} o.velocity translator speed (instantaneous or amplitude), m/s
 * @param {number} o.Rg phase resistance, Ω @param {number} o.RL load resistance per phase, Ω
 * @returns {{emf:number, current:number, force:number, damping:number, outputPower:number, copperLoss:number,
 *            efficiency:number, meanOutputSinusoidal:number}} V, A (peak), N, N·s/m, W, W, –, W
 */
export function linearGenerator({ ke, velocity, Rg, RL }) {
  pos(ke, 'ke'); fin(velocity, 'velocity'); nonNeg(Rg, 'Rg'); pos(RL, 'RL');
  const emf = ke * velocity;
  const current = emf / (Rg + RL);
  const damping = (1.5 * ke * ke) / (Rg + RL);
  const outputPower = 1.5 * current * current * RL;
  return { emf, current, force: damping * velocity, damping, outputPower, copperLoss: 1.5 * current * current * Rg,
    efficiency: RL / (Rg + RL), meanOutputSinusoidal: 0.5 * outputPower };
}

/**
 * Electrical side of a hydraulic power take-off as an efficiency chain with a part-load characteristic:
 * Pout = Pin·ηpump·ηmotor(x)·ηgenerator(x)·ηconverter, x = Pin/Prated, where the hydraulic motor and the
 * generator have η(x) = x/(x + l0 + l1·x²): a constant loss l0 and a load-squared loss l1 (fractions of rated
 * output) with l0 + l1 = 1/ηrated − 1, so that η(1) = ηrated. The split between l0 and l1 is equipment data
 * (default half each). A simple efficiency model, not a hydraulic circuit simulation.
 * @param {object} o
 * @param {number} o.power absorbed hydraulic input power, W @param {number} o.ratedPower W
 * @param {number} [o.pumpEfficiency=0.95] @param {number} [o.motorEfficiency=0.92] @param {number} [o.generatorEfficiency=0.95]
 * @param {number} [o.converterEfficiency=0.97] @param {number} [o.constantLossShare=0.5] share of the rated loss of
 *        motor and generator that is independent of load
 * @returns {{outputPower:number, efficiency:number, stages:{pump:number, motor:number, generator:number, converter:number}}}
 */
export function hydraulicPto({ power, ratedPower, pumpEfficiency = 0.95, motorEfficiency = 0.92, generatorEfficiency = 0.95,
  converterEfficiency = 0.97, constantLossShare = 0.5 }) {
  nonNeg(power, 'power'); pos(ratedPower, 'ratedPower');
  for (const [v, n] of [[pumpEfficiency, 'pumpEfficiency'], [motorEfficiency, 'motorEfficiency'], [generatorEfficiency, 'generatorEfficiency'], [converterEfficiency, 'converterEfficiency']]) {
    if (!(v > 0 && v <= 1)) throw new RangeError(`${n}: must be in (0, 1]`);
  }
  if (!(constantLossShare >= 0 && constantLossShare <= 1)) throw new RangeError('constantLossShare: must be in [0, 1]');
  const partLoad = (etaRated, x) => {
    if (x <= 0) return 0;
    const loss = 1 / etaRated - 1; // rated loss per unit rated output
    const l0 = constantLossShare * loss; const l1 = (1 - constantLossShare) * loss;
    return x / (x + l0 + l1 * x * x);
  };
  const x = power / ratedPower;
  const stages = { pump: pumpEfficiency, motor: partLoad(motorEfficiency, x), generator: partLoad(generatorEfficiency, x), converter: converterEfficiency };
  const efficiency = stages.pump * stages.motor * stages.generator * stages.converter;
  return { outputPower: power * efficiency, efficiency, stages };
}
