/**
 * Incremental nonlinear static (pushover) collapse analysis of frames with concentrated plasticity, intended
 * for ultimate-strength screening of tubular space frames (jackets, towers, module frames).
 *
 * Method — event-to-event plastic-hinge analysis (see e.g. Chen & Sohal, Plastic Design and Second-Order Analysis
 * of Steel Frames, 1995):
 *  - A load pattern is scaled by the load factor λ; between events the response is linear with the current
 *    tangent stiffness (plus geometric stiffness of the current member axial forces when `pDelta` is on).
 *  - Plastic hinges form at member ends when the cross-section interaction reaches unity. Tubes use
 *    cos(π/2 · N/Np) − M/Mp = 0 with M = √(My² + Mz²) (fully plastic stress distribution of a thin-walled circular
 *    tube). Other sections use the linear lower bound |N|/Np + |My|/Mpy + |Mz|/Mpz = 1.
 *    A formed hinge is modelled as a bending release (both axes) that keeps its moment — elastic-perfectly-plastic,
 *    no hardening, no unloading, and the hinge moment is not moved along the yield surface when N changes later.
 *  - Tension members yield at Np (or at a joint capacity, if lower) and then elongate at constant force.
 *  - Compression members fail at the column-buckling capacity Pc (ISO 19902 §13.2.3 curve by default). The
 *    post-buckling branch is a documented simple model: the member force drops to `postBucklingRatio`·Pc and stays
 *    constant; the released force is redistributed at constant λ through an inner event-to-event stage, so a
 *    cascade of failures (dynamic collapse is not modelled) is followed.
 *  - Collapse is declared when the tangent stiffness becomes singular or indefinite (mechanism / instability) or a
 *    translational degree of freedom is left without stiffness.
 * Results: load factor vs control displacement, base shear, failure sequence, ultimate capacity, reserve strength
 * ratio RSR = ultimate / design load, damaged strength and redundancy measures, member-removal screening.
 * SI units.
 */

import {
  buildFrameSystem, factorFrame, solveFrame, resolveLoads, loadParts,
} from './frame3d.js';
import { frameGeometricStiffness } from './framebuckling.js';
import { isoTubularStrengths } from './tubular.js';

/**
 * Plastic capacities of a member used by the pushover analysis. `member.capacity` ({Np, Mpy, Mpz, Pc, Pt, Pj})
 * overrides any of the defaults: Np = A fy; Mp = Wpl fy; tension limit Pt = Np; compression limit Pc from the
 * ISO 19902 column curve fc = (1 − 0.278 λ²) fyc (λ ≤ 1.34), 0.9 fyc/λ² otherwise, with λ = (K L/(π r))√(fyc/E),
 * K = member.K (default 1) and L = member.bucklingLength (default the member length); local buckling reduction
 * fyc for tubes per §13.2.3.3. Pj is an optional joint capacity that caps both Pt and Pc.
 * @param {object} elem element record from buildFrameSystem
 * @returns {{Np:number, Mpy:number, Mpz:number, Pt:number, Pc:number, Pj:number, tube:boolean}} forces [N], moments [N·m]
 */
export function memberPlasticCapacity(elem) {
  const sec = elem.sec; const mat = elem.mat; const mem = elem.member;
  const fy = mat.fy;
  const cap = mem.capacity || {};
  const tube = sec.type === 'tube';
  const Np = cap.Np !== undefined ? cap.Np : sec.A * fy;
  const wy = sec.Wply || sec.Wely || (sec.cz ? sec.Iy / sec.cz : NaN);
  const wz = sec.Wplz || sec.Welz || (sec.cy ? sec.Iz / sec.cy : NaN);
  const Mpy = cap.Mpy !== undefined ? cap.Mpy : (cap.Mp !== undefined ? cap.Mp : wy * fy);
  const Mpz = cap.Mpz !== undefined ? cap.Mpz : (cap.Mp !== undefined ? cap.Mp : wz * fy);
  let Pc = cap.Pc;
  if (Pc === undefined) {
    const K = mem.K === undefined ? 1 : mem.K;
    const Lb = mem.bucklingLength === undefined ? elem.L : mem.bucklingLength;
    if (tube) Pc = isoTubularStrengths({ D: sec.D, t: sec.t, fy, E: mat.E, L: Lb, K }).Pc;
    else {
      const r = Math.sqrt(Math.min(sec.Iy, sec.Iz) / sec.A);
      const lam = K * Lb / (Math.PI * r) * Math.sqrt(fy / mat.E);
      Pc = sec.A * fy * (lam <= 1.34 ? 1 - 0.278 * lam * lam : 0.9 / (lam * lam));
    }
  }
  const Pj = cap.Pj !== undefined ? cap.Pj : (mem.jointCapacity !== undefined ? mem.jointCapacity : Infinity);
  const Pt = cap.Pt !== undefined ? cap.Pt : Np;
  if (![Np, Mpy, Mpz, Pc, Pt].every((v) => v > 0)) throw new RangeError(`member ${elem.index}: plastic capacities must be positive`);
  return { Np, Mpy, Mpz, Pt, Pc, Pj, tube };
}

/**
 * Plastic interaction function of a member cross-section: positive inside the yield surface, zero on it.
 * Tubes: g = cos(π/2 · N/Np) − √(My² + Mz²)/Mp. Other sections: g = 1 − |N|/Np − |My|/Mpy − |Mz|/Mpz.
 * @param {{Np:number, Mpy:number, Mpz:number, tube:boolean}} cap capacities [N, N·m]
 * @param {number} N axial force [N]
 * @param {number} My bending moment [N·m]
 * @param {number} Mz bending moment [N·m]
 * @returns {number} g [-]
 */
export function plasticInteraction(cap, N, My, Mz) {
  const n = Math.min(Math.abs(N) / cap.Np, 2);
  if (cap.tube) return Math.cos(Math.PI / 2 * n) - Math.hypot(My, Mz) / cap.Mpy;
  return 1 - n - Math.abs(My) / cap.Mpy - Math.abs(Mz) / cap.Mpz;
}

function withReleases(mem, hingeA, hingeB, axial) {
  if (!hingeA && !hingeB && !axial) return mem;
  const flags = (r) => {
    const out = [0, 0, 0, 0, 0, 0];
    if (!r) return out;
    if (r === 'pin' || r === 'pinned') { out[4] = 1; out[5] = 1; return out; }
    if (typeof r === 'string') {
      const map = { ux: 0, uy: 1, uz: 2, rx: 3, ry: 4, rz: 5, n: 0, vy: 1, vz: 2, t: 3, my: 4, mz: 5 };
      for (const tok of r.split(/[\s,]+/)) if (tok) out[map[tok.toLowerCase()]] = 1;
      return out;
    }
    for (let i = 0; i < 6; i++) out[i] = r[i] ? 1 : 0;
    return out;
  };
  const a = flags(mem.releases && mem.releases.a);
  const b = flags(mem.releases && mem.releases.b);
  if (hingeA) { a[4] = 1; a[5] = 1; }
  if (hingeB) { b[4] = 1; b[5] = 1; }
  if (axial && !a[0]) b[0] = 1;
  return { ...mem, releases: { a, b } };
}

/**
 * Event-to-event pushover analysis.
 * @param {object} model frame model with load cases (frame3d format). Optional member fields: capacity {Np, Mp|Mpy, Mpz, Pc, Pt, Pj},
 *   K, bucklingLength, jointCapacity.
 * @param {object} o
 * @param {string} o.pattern load case / combination scaled by the load factor λ
 * @param {string} [o.base] load case / combination applied first at factor 1 and held (permanent loads)
 * @param {{node:number, dof:(number|string)}} [o.control] displacement reported in the curve; default: the translational dof with
 *   the largest elastic response to the pattern
 * @param {boolean} [o.pDelta=true] include the geometric stiffness of the current axial forces
 * @param {number} [o.postBucklingRatio=0.4] residual compression capacity after buckling as a fraction of Pc [-], 0..1
 * @param {number} [o.designFactor=1] load factor of the design load level (RSR = λult/designFactor)
 * @param {number} [o.maxEvents=500]
 * @param {number} [o.maxLoadFactor=Infinity] stop when λ reaches this value
 * @param {number} [o.maxStep] largest λ increment per step (default: unlimited without P-Δ; one fifth of the first event with P-Δ)
 * @param {boolean} [o.shearDeformation=true]
 * @returns {{curve:{loadFactor:number[], displacement:number[], baseShear:number[]}, events:Array<{loadFactor:number, member:number,
 *   type:string, displacement:number, stage:string}>, ultimateLoadFactor:number, firstEventLoadFactor:number, RSR:number, status:string,
 *   baseShearRef:number, ultimateBaseShear:number, redundancyFactor:number, memberStates:Array<{axial:string, hingeA:boolean, hingeB:boolean,
 *   N:number}>, displacements:Float64Array, capacities:object[], control:{node:number, dof:number}}}
 *   status: 'collapse' | 'collapse-under-base-load' | 'max-load-factor' | 'max-events' | 'no-further-events';
 *   redundancyFactor = λult / λ at the first member failure
 */
export function pushoverAnalysis(model, {
  pattern, base, control, pDelta = true, postBucklingRatio = 0.4, designFactor = 1, maxEvents = 500, maxLoadFactor = Infinity,
  maxStep, shearDeformation = true,
} = {}) {
  if (typeof pattern !== 'string') throw new RangeError('pattern must name a load case or combination');
  if (!(postBucklingRatio >= 0 && postBucklingRatio <= 1)) throw new RangeError('postBucklingRatio must lie in [0, 1]');
  if (!(designFactor > 0)) throw new RangeError('designFactor must be positive');
  const sys0 = buildFrameSystem(model, { shearDeformation });
  const nEl = sys0.elems.length;
  const caps = sys0.elems.map(memberPlasticCapacity);
  const st = sys0.elems.map((el) => ({
    axial: 'elastic', hingeA: false, hingeB: false, f: new Float64Array(12),
    canHingeA: !el.truss && !(el.relIdx.includes(4) && el.relIdx.includes(5)),
    canHingeB: !el.truss && !(el.relIdx.includes(10) && el.relIdx.includes(11)),
  }));
  const u = new Float64Array(sys0.nDof);
  const patternParts = loadParts(model, pattern);
  const DOFS = { ux: 0, uy: 1, uz: 2, rx: 3, ry: 4, rz: 5 };
  let ctrl = null;
  if (control) {
    const k = typeof control.dof === 'string' ? DOFS[control.dof] : control.dof;
    if (!Number.isInteger(control.node) || control.node < 0 || control.node >= sys0.nNode || !(k >= 0 && k < 6)) throw new RangeError('control needs a valid node and dof');
    ctrl = 6 * control.node + k;
  }
  const axialOf = (f) => 0.5 * (-f[0] + f[6]);

  // tangent system for the current member states
  const tangent = () => {
    const members = model.members.map((m, e) => withReleases(m, st[e].hingeA, st[e].hingeB, st[e].axial !== 'elastic'));
    const sys = buildFrameSystem({ ...model, members }, { shearDeformation });
    for (const d of sys.autoRestrained) {
      if (d % 6 < 3 && !sys0.auto[d]) {
        const err = new Error(`Translational dof ${d} has lost all stiffness.`);
        err.mechanism = true;
        throw err;
      }
    }
    let kg = null;
    if (pDelta) kg = frameGeometricStiffness(sys, st.map((s) => axialOf(s.f)));
    const fac = factorFrame(sys, kg ? { add: kg.global } : {});
    if (fac.negativePivots > 0) {
      const err = new Error('Tangent stiffness is indefinite (global instability).');
      err.mechanism = true;
      throw err;
    }
    return { sys, kg, fac };
  };
  const unit = (tg, loads) => solveFrame(tg.sys, tg.fac, loads, tg.kg ? { add: tg.kg.global, addLocal: tg.kg.local } : {});

  // smallest factor at which any member reaches a limit along the direction of the unit response
  const nextEvent = (res) => {
    let dmin = Infinity;
    const cand = [];
    const consider = (d, e, type) => {
      if (!(d >= 0)) return;
      cand.push({ d, e, type });
      if (d < dmin) dmin = d;
    };
    for (let e = 0; e < nEl; e++) {
      const s = st[e]; const cap = caps[e];
      const f0 = s.f; const f1 = res.memberForces[e].f;
      if (s.axial === 'elastic') {
        const N0 = axialOf(f0); const n1 = axialOf(f1);
        const Pt = Math.min(cap.Pt, cap.Pj); const Pc = Math.min(cap.Pc, cap.Pj);
        const scale = 1e-12 * Math.max(Pt, Pc);
        if (n1 > scale * 1e-3) consider(Math.max(0, (Pt - N0) / n1), e, cap.Pj < cap.Pt ? 'joint' : 'tension-yield');
        else if (n1 < -scale * 1e-3) consider(Math.max(0, (N0 + Pc) / -n1), e, cap.Pj < cap.Pc ? 'joint' : 'buckling');
      }
      for (const end of [0, 1]) {
        if (end === 0 ? (s.hingeA || !s.canHingeA) : (s.hingeB || !s.canHingeB)) continue;
        const o = end === 0 ? 0 : 6; const sg = end === 0 ? -1 : 1;
        const N0 = sg * f0[o]; const My0 = sg * f0[o + 4]; const Mz0 = sg * f0[o + 5];
        const n1 = sg * f1[o]; const my1 = sg * f1[o + 4]; const mz1 = sg * f1[o + 5];
        const mScale = Math.max(Math.abs(my1), Math.abs(mz1)) / Math.min(cap.Mpy, cap.Mpz) + Math.abs(n1) / cap.Np;
        if (!(mScale > 0)) continue;
        const g = (d) => plasticInteraction(cap, N0 + d * n1, My0 + d * my1, Mz0 + d * mz1);
        const g0 = g(0);
        // the interaction function is concave along the load path: at most one outward crossing
        let hi = 1e-9 / mScale;
        if (g0 <= 1e-10) {
          if (g(hi) < g0) consider(0, e, end === 0 ? 'hinge-A' : 'hinge-B');
          continue;
        }
        let found = false;
        const cap2 = 4 / mScale; // beyond this the section forces alone exceed twice the capacity
        while (hi <= cap2) { if (g(hi) < 0) { found = true; break; } hi *= 2; }
        if (!found) { if (g(cap2) < 0) { hi = cap2; found = true; } else continue; }
        let lo = hi / 2 > 1e-9 / mScale ? hi / 2 : 0;
        if (g(lo) < 0) lo = 0;
        for (let it = 0; it < 80; it++) { const mid = 0.5 * (lo + hi); if (g(mid) > 0) lo = mid; else hi = mid; }
        consider(hi, e, end === 0 ? 'hinge-A' : 'hinge-B');
      }
    }
    const tol = dmin * 1e-9 + 1e-14;
    return { d: dmin, hits: cand.filter((c) => c.d <= dmin + tol) };
  };

  const curve = { loadFactor: [0], displacement: [0], baseShear: [0] };
  const events = [];
  let lambda = 0; let status = null; let baseShearRef = 0; let lambdaMax = 0;
  const pending = []; // redistribution loads after buckling: { F, remaining }
  const zeroItems = () => sys0.elems.map(() => []);
  let stageName = base ? 'base' : 'pattern';
  let baseRemaining = base ? 1 : 0;
  let stepCap = maxStep;
  let nEvents = 0;
  let guard = 0;

  const apply = (res, d) => {
    for (let g = 0; g < sys0.nDof; g++) u[g] += d * res.displacements[g];
    for (let e = 0; e < nEl; e++) { const f = st[e].f; const f1 = res.memberForces[e].f; for (let i = 0; i < 12; i++) f[i] += d * f1[i]; }
  };
  const record = () => {
    curve.loadFactor.push(lambda);
    curve.displacement.push(ctrl === null ? 0 : u[ctrl]);
    curve.baseShear.push(lambda * baseShearRef);
  };
  const trigger = (hits, stage) => {
    for (const h of hits) {
      const s = st[h.e]; const cap = caps[h.e]; const el = sys0.elems[h.e];
      if (h.type === 'hinge-A') { if (s.hingeA) continue; s.hingeA = true; }
      else if (h.type === 'hinge-B') { if (s.hingeB) continue; s.hingeB = true; }
      else {
        if (s.axial !== 'elastic') continue;
        if (h.type === 'buckling') {
          s.axial = 'buckled';
          const drop = (1 - postBucklingRatio) * Math.min(cap.Pc, cap.Pj);
          if (drop > 0) {
            // the member pushes its end nodes apart with Pc; after buckling only with ρ·Pc
            const F = new Float64Array(sys0.nDof);
            const ex = el.R[0];
            for (let k = 0; k < 3; k++) { F[6 * el.a + k] += drop * ex[k]; F[6 * el.b + k] -= drop * ex[k]; }
            pending.push({ F, remaining: 1 });
            const res = -postBucklingRatio * Math.min(cap.Pc, cap.Pj);
            s.f[0] = -res; s.f[6] = res;
          }
        } else s.axial = h.type === 'joint' ? 'joint' : 'yielded';
      }
      nEvents++;
      events.push({ loadFactor: lambda, member: h.e, type: h.type, displacement: ctrl === null ? 0 : u[ctrl], stage });
    }
  };

  while (!status) {
    if (++guard > 20 * maxEvents + 1000) { status = 'max-events'; break; }
    let tg;
    try {
      tg = tangent();
    } catch (err) {
      if (!err.mechanism) throw err;
      status = stageName === 'base' ? 'collapse-under-base-load' : 'collapse';
      break;
    }
    let loads; let limit; let mode;
    if (pending.length) {
      loads = { F: pending[0].F, items: zeroItems(), prescribed: new Float64Array(sys0.nDof) };
      limit = pending[0].remaining; mode = 'redistribution';
    } else if (baseRemaining > 0) {
      loads = resolveLoads(tg.sys, loadParts(model, base)); limit = baseRemaining; mode = 'base';
    } else {
      loads = resolveLoads(tg.sys, patternParts); limit = maxLoadFactor - lambda; mode = 'pattern';
      stageName = 'pattern';
    }
    let res;
    try {
      res = unit(tg, loads);
    } catch (err) {
      if (!err.mechanism) throw err;
      status = stageName === 'base' ? 'collapse-under-base-load' : 'collapse';
      break;
    }
    if (mode === 'pattern' && baseShearRef === 0) {
      const a = res.appliedSum;
      const h = Math.hypot(a[0], a[1]);
      baseShearRef = h > 1e-9 * Math.abs(a[2]) ? h : Math.abs(a[2]);
      if (ctrl === null) {
        let best = 0;
        for (let g = 0; g < sys0.nDof; g++) if (g % 6 < 3 && Math.abs(res.displacements[g]) > best) { best = Math.abs(res.displacements[g]); ctrl = g; }
      }
    }
    const ev = nextEvent(res);
    let d = ev.d;
    let hits = ev.hits;
    if (mode === 'pattern' && pDelta && stepCap === undefined && Number.isFinite(d) && d > 0) stepCap = d / 5;
    let capped = false;
    if (mode === 'pattern' && stepCap !== undefined && d > stepCap) { d = stepCap; capped = true; }
    if (d > limit) { d = limit; capped = true; }
    if (!Number.isFinite(d)) { status = 'no-further-events'; break; }
    apply(res, d);
    if (mode === 'pattern') { lambda += d; if (lambda > lambdaMax) lambdaMax = lambda; record(); }
    else if (mode === 'base') { baseRemaining -= d; if (baseRemaining < 1e-12) baseRemaining = 0; }
    else { pending[0].remaining -= d; if (pending[0].remaining < 1e-12) pending.shift(); record(); }
    if (!capped) trigger(hits, mode);
    else if (mode === 'pattern' && lambda >= maxLoadFactor * (1 - 1e-12)) status = 'max-load-factor';
    if (nEvents >= maxEvents) status = 'max-events';
  }
  const firstEvent = events.find((e) => e.stage !== 'base');
  const memberStates = st.map((s) => ({ axial: s.axial, hingeA: s.hingeA, hingeB: s.hingeB, N: axialOf(s.f) }));
  return {
    curve, events, ultimateLoadFactor: lambdaMax, firstEventLoadFactor: firstEvent ? firstEvent.loadFactor : lambdaMax,
    RSR: lambdaMax / designFactor, status, baseShearRef, ultimateBaseShear: lambdaMax * baseShearRef,
    redundancyFactor: firstEvent && firstEvent.loadFactor > 0 ? lambdaMax / firstEvent.loadFactor : 1,
    memberStates, displacements: u, capacities: caps, control: ctrl === null ? null : { node: Math.floor(ctrl / 6), dof: ctrl % 6 },
  };
}

/**
 * Copy of a model with members removed; member loads on removed members are dropped and the member indices of
 * the remaining member loads are renumbered.
 * @param {object} model frame model
 * @param {number[]} remove member indices to delete
 * @returns {{model:object, indexMap:Int32Array}} new model and old → new member index map (−1 for removed members)
 */
export function removeMembers(model, remove) {
  const drop = new Set(remove);
  for (const r of drop) if (!Number.isInteger(r) || r < 0 || r >= model.members.length) throw new RangeError(`member index ${r} is out of range`);
  const indexMap = new Int32Array(model.members.length).fill(-1);
  const members = [];
  model.members.forEach((m, i) => { if (!drop.has(i)) { indexMap[i] = members.length; members.push(m); } });
  if (members.length === 0) throw new RangeError('cannot remove every member');
  const out = { ...model, members };
  if (model.loadCases) {
    out.loadCases = model.loadCases.map((lc) => {
      if (!lc.member) return lc;
      return { ...lc, member: lc.member.filter((ml) => indexMap[ml.member] >= 0).map((ml) => ({ ...ml, member: indexMap[ml.member] })) };
    });
  }
  return { model: out, indexMap };
}

/**
 * Residual (damaged) strength with named members removed, compared with the intact structure.
 * Damaged strength ratio DSR = λult,damaged/λult,intact (residual resistance factor); RSR of the damaged
 * structure = λult,damaged/designFactor; redundancy factor of the intact structure = λult / λ at first member failure.
 * @param {object} model frame model
 * @param {object} o pushoverAnalysis options plus
 * @param {number[]} o.remove indices of the members removed (damaged / severed)
 * @returns {{intact:object, damaged:object, damagedStrengthRatio:number, RSRintact:number, RSRdamaged:number, redundancyFactor:number}}
 */
export function damagedStrength(model, { remove, ...opts } = {}) {
  if (!Array.isArray(remove) || remove.length === 0) throw new RangeError('remove must list at least one member');
  const intact = pushoverAnalysis(model, opts);
  const dm = removeMembers(model, remove).model;
  let damaged;
  try {
    damaged = pushoverAnalysis(dm, opts);
  } catch (err) {
    if (!err.mechanism) throw err;
    damaged = { ultimateLoadFactor: 0, RSR: 0, status: 'collapse', events: [], curve: { loadFactor: [0], displacement: [0], baseShear: [0] } };
  }
  return {
    intact, damaged, damagedStrengthRatio: intact.ultimateLoadFactor > 0 ? damaged.ultimateLoadFactor / intact.ultimateLoadFactor : 0,
    RSRintact: intact.RSR, RSRdamaged: damaged.RSR, redundancyFactor: intact.redundancyFactor,
  };
}

/**
 * Progressive-collapse screening by notional member removal (alternate load path method, cf. GSA 2003 /
 * UFC 4-023-03): each listed member is removed in turn and the remaining structure is pushed with the given load
 * pattern; it survives when its ultimate load factor is at least the dynamic amplification factor applied to the
 * loads (2.0 is the classical upper bound for sudden removal in a static analysis; pass 1.0 for a static check).
 * @param {object} model frame model
 * @param {object} o pushoverAnalysis options (pattern = the load combination present at removal) plus
 * @param {number[]} [o.members] members to remove one at a time (default: all)
 * @param {number} [o.daf=2.0] dynamic amplification factor on the loads [-]
 * @returns {{cases:Array<{member:number, ultimateLoadFactor:number, demand:number, uc:number, survives:boolean, status:string, failures:number}>,
 *   critical:object, allSurvive:boolean}} uc = daf/λult
 */
export function progressiveCollapse(model, { members, daf = 2.0, ...opts } = {}) {
  if (!(daf >= 1)) throw new RangeError('daf must be at least 1');
  const list = members || model.members.map((m, i) => i);
  const cases = [];
  for (const m of list) {
    const dm = removeMembers(model, [m]).model;
    let r;
    try {
      r = pushoverAnalysis(dm, { ...opts, maxLoadFactor: opts.maxLoadFactor === undefined ? 4 * daf : opts.maxLoadFactor });
    } catch (err) {
      if (!err.mechanism) throw err;
      r = { ultimateLoadFactor: 0, status: 'collapse', events: [] };
    }
    const lam = r.ultimateLoadFactor;
    cases.push({
      member: m, ultimateLoadFactor: lam, demand: daf, uc: lam > 0 ? daf / lam : Infinity, survives: lam >= daf, status: r.status,
      failures: r.events.length,
    });
  }
  let critical = cases[0];
  for (const c of cases) if (c.uc > critical.uc) critical = c;
  return { cases, critical, allSurvive: cases.every((c) => c.survives) };
}
