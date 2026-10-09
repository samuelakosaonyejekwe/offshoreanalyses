// Floating-solar platform engineering: buoyancy and freeboard, island stability, wind/wave/current loads,
// mooring load sharing and water-level slack, hinged-raft kinematics and connector screening, evaporation and
// cooling benefits, floating-cable slack and site-exposure screening.
//
// Units: SI (N, m, kg, Pa, rad); angles in degrees where the name ends in Deg.

const G = 9.80665;
const DEG = Math.PI / 180;

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function finitePos(v, name) {
  need(Number.isFinite(v) && v > 0, `${name} must be a positive finite number`);
}

/**
 * Buoyancy, draft and freeboard of a float module (pontoon or HDPE float treated as a prism of height h and
 * effective waterplane area V/h): draft = m / (ρ A_wp), freeboard = h − draft, reserve buoyancy = ρ V/m − 1.
 * Three load cases are reported: dead load, dead + maintenance live load, dead + live + snow.
 * Archimedes' principle; acceptance limits (minimum freeboard, reserve-buoyancy ratio) are project criteria —
 * DNV-RP-0584 (2021) gives design guidance for floating PV.
 * @param {object} o
 * @param {number} o.floatVolume enclosed volume of the float(s) supporting the tributary load (m³)
 * @param {number} o.floatHeight float height h (m)
 * @param {number} o.deadMass supported permanent mass incl. float, module, frame, cables (kg)
 * @param {number} [o.liveMass=0] maintenance live load (kg)
 * @param {number} [o.snowMass=0] snow load on the tributary area (kg)
 * @param {number} [o.rho=1000] water density (kg/m³)
 * @param {number} [o.minFreeboard=0] required freeboard (m)
 * @returns {{cases:{name:string, mass:number, draft:number, freeboard:number, reserveBuoyancy:number,
 *   utilisation:number, ok:boolean}[], waterplaneArea:number}} utilisation = weight / total buoyancy (-)
 */
export function floatBuoyancy({ floatVolume, floatHeight, deadMass, liveMass = 0, snowMass = 0, rho = 1000, minFreeboard = 0 }) {
  finitePos(floatVolume, 'floatVolume');
  finitePos(floatHeight, 'floatHeight');
  finitePos(deadMass, 'deadMass');
  finitePos(rho, 'rho');
  need(liveMass >= 0 && snowMass >= 0, 'liveMass and snowMass must be >= 0');
  const Awp = floatVolume / floatHeight;
  const mk = (name, mass) => {
    const draft = mass / (rho * Awp);
    const freeboard = floatHeight - draft;
    return { name, mass, draft, freeboard, reserveBuoyancy: (rho * floatVolume) / mass - 1, utilisation: mass / (rho * floatVolume), ok: freeboard >= minFreeboard && freeboard > 0 };
  };
  return {
    cases: [mk('dead', deadMass), mk('dead+live', deadMass + liveMass), mk('dead+live+snow', deadMass + liveMass + snowMass)],
    waterplaneArea: Awp,
  };
}

/**
 * Intact (initial) transverse stability of a rigid float-island segment:
 *   I_T = Σ (a_i y_i² + i_i)  or  φ_wp L B³/12 for uniformly distributed floats with waterplane fill ratio φ_wp,
 *   BM = I_T / ∇, ∇ = m/ρ,  KB = draft/2,  GM = KB + BM − KG,
 *   heel under a moment M: sin φ = M / (m g GM).
 * Heeling moment = wind force × lever + person weight × offset. Small-angle hydrostatics (PNA vol. 1); valid
 * while no float edge submerges or emerges.
 * @param {object} o
 * @param {number} o.mass segment mass (kg)
 * @param {number} o.kg height of the centre of gravity above the float bottom KG (m)
 * @param {number} o.draft float draft (m)
 * @param {number} [o.length] segment length L (m), with width and fillRatio
 * @param {number} [o.width] segment width B (m), heeling about the longitudinal axis
 * @param {number} [o.fillRatio=1] waterplane area / (L B) (-)
 * @param {{y:number, area:number, ownInertia?:number}[]} [o.floats] explicit floats: offset from the centreline
 *   (m), waterplane area (m²), own waterplane inertia (m⁴); overrides length/width/fillRatio
 * @param {number} [o.windForce=0] horizontal wind force (N)
 * @param {number} [o.windLever=0] height of the wind force above the mooring/lateral-resistance level (m)
 * @param {number} [o.personMass=0] maintenance personnel mass (kg)
 * @param {number} [o.personOffset=0] transverse offset of personnel from the centreline (m)
 * @param {number} [o.rho=1000]
 * @returns {{waterplaneInertia:number, BM:number, KB:number, GM:number, heelingMoment:number, heelDeg:number,
 *   stable:boolean}}
 */
export function islandStability({ mass, kg, draft, length, width, fillRatio = 1, floats, windForce = 0, windLever = 0, personMass = 0, personOffset = 0, rho = 1000 }) {
  finitePos(mass, 'mass');
  finitePos(draft, 'draft');
  need(Number.isFinite(kg), 'kg must be finite');
  let I;
  if (floats) {
    need(Array.isArray(floats) && floats.length > 0, 'floats must be a non-empty array');
    I = 0;
    for (const f of floats) I += f.area * f.y * f.y + (f.ownInertia ?? 0);
  } else {
    finitePos(length, 'length');
    finitePos(width, 'width');
    need(fillRatio > 0 && fillRatio <= 1, 'fillRatio must be within (0, 1]');
    I = (fillRatio * length * width ** 3) / 12;
  }
  const BM = I / (mass / rho);
  const KB = draft / 2;
  const GM = KB + BM - kg;
  const M = windForce * windLever + personMass * G * personOffset;
  const s = GM > 0 ? M / (mass * G * GM) : NaN;
  return {
    waterplaneInertia: I,
    BM,
    KB,
    GM,
    heelingMoment: M,
    heelDeg: GM > 0 && Math.abs(s) <= 1 ? Math.asin(s) / DEG : NaN,
    stable: GM > 0,
  };
}

/**
 * Wind load on rows of tilted PV panels. The net pressure on a panel acts normal to its surface,
 * F_N = q c_N A with q = ½ ρ V²; resolved into drag F_N sin β and lift F_N cos β (uplift for wind on the rear
 * face, downforce for wind on the front). The windward row takes the full load; each inner row is multiplied by
 * a shelter factor.
 * The net normal-force coefficient c_N and the shelter factor are project inputs without defaults: take them
 * from the governing code for monoslope/free-standing canopies (EN 1991-1-4 §7.3, monopitch canopies; ASCE 7 provisions for
 * open monoslope roofs and solar arrays) or from wind-tunnel tests of the float system.
 * @param {object} o
 * @param {number} o.windSpeed design wind speed at panel height (m/s)
 * @param {number} o.panelArea area of one row's panels (m²)
 * @param {number} o.tiltDeg panel tilt β (deg)
 * @param {number} o.cN net normal-force coefficient (-)
 * @param {number} [o.rows=1] number of rows along the wind
 * @param {number} [o.shelterFactor=1] load factor of inner rows relative to the first (-), 0–1
 * @param {number} [o.rho=1.225] air density (kg/m³)
 * @returns {{dynamicPressure:number, normalForcePerRow:number, drag:number, lift:number, dragFirstRow:number,
 *   liftFirstRow:number}} Pa and N; drag and lift are totals over all rows (lift magnitude)
 */
export function panelWindLoad({ windSpeed, panelArea, tiltDeg, cN, rows = 1, shelterFactor = 1, rho = 1.225 }) {
  need(windSpeed >= 0, 'windSpeed must be >= 0');
  finitePos(panelArea, 'panelArea');
  need(Number.isFinite(cN), 'cN must be given (net normal-force coefficient from the governing wind code)');
  need(tiltDeg >= 0 && tiltDeg <= 90, 'tiltDeg must be within [0, 90]');
  need(Number.isInteger(rows) && rows >= 1, 'rows must be an integer >= 1');
  need(shelterFactor >= 0 && shelterFactor <= 1, 'shelterFactor must be within [0, 1]');
  const q = 0.5 * rho * windSpeed * windSpeed;
  const FN = q * Math.abs(cN) * panelArea;
  const mult = 1 + (rows - 1) * shelterFactor;
  const d1 = FN * Math.sin(tiltDeg * DEG);
  const l1 = FN * Math.cos(tiltDeg * DEG);
  return { dynamicPressure: q, normalForcePerRow: FN, drag: d1 * mult, lift: l1 * mult, dragFirstRow: d1, liftFirstRow: l1 };
}

function waveNumber(omega, h) {
  const k0 = (omega * omega) / G;
  if (!Number.isFinite(h) || k0 * h > 40) return k0;
  let k = k0 / Math.sqrt(Math.tanh(k0 * h));
  for (let i = 0; i < 60; i++) {
    const th = Math.tanh(k * h);
    const dk = (G * k * th - omega * omega) / (G * th + G * k * h * (1 - th * th));
    k -= dk;
    if (Math.abs(dk) < 1e-14 * k) break;
  }
  return k;
}

/**
 * Horizontal wave force on a small float in the long-wave (Morison / Froude–Krylov) regime: inertia amplitude
 * F_I = ρ C_m V a_x and drag amplitude F_D = ½ ρ C_d A u², with surface particle kinematics of a linear wave,
 * u = ω (H/2)/tanh(kh), a_x = ω u. Maximum of the combined force over a cycle: F_I if F_I ≥ 2 F_D, otherwise
 * F_D + F_I²/(4 F_D). Valid for float dimension ≪ wavelength (D/λ < 0.2); reference DNV-RP-C205 §6.
 * @param {object} o
 * @param {number} o.H wave height (m)
 * @param {number} o.T wave period (s)
 * @param {number} [o.depth=Infinity] water depth (m)
 * @param {number} o.submergedVolume displaced volume of the float (m³)
 * @param {number} o.projectedArea submerged area normal to the wave direction (m²)
 * @param {number} [o.Cm=2] inertia coefficient 1 + C_a (-)
 * @param {number} [o.Cd=1] drag coefficient (-)
 * @param {number} [o.rho=1000]
 * @returns {{inertia:number, drag:number, max:number, velocity:number, acceleration:number, wavelength:number}}
 */
export function floatWaveLoad({ H, T, depth = Infinity, submergedVolume, projectedArea, Cm = 2, Cd = 1, rho = 1000 }) {
  need(H >= 0, 'H must be >= 0');
  finitePos(T, 'T');
  need(submergedVolume >= 0 && projectedArea >= 0, 'submergedVolume and projectedArea must be >= 0');
  const w = (2 * Math.PI) / T;
  const k = waveNumber(w, depth);
  const th = Number.isFinite(depth) ? Math.tanh(k * depth) : 1;
  const u = (w * H) / 2 / th;
  const a = w * u;
  const FI = rho * Cm * submergedVolume * a;
  const FD = 0.5 * rho * Cd * projectedArea * u * u;
  const max = FD === 0 || FI >= 2 * FD ? FI : FD + (FI * FI) / (4 * FD);
  return { inertia: FI, drag: FD, max, velocity: u, acceleration: a, wavelength: (2 * Math.PI) / k };
}

/**
 * Mean wave-drift force on a float island from wave reflection (Maruo 1960, deep water, long-crested waves):
 * regular wave F̄ = ½ ρ g a² R² per unit width; irregular sea F̄ = ρ g R² Hs²/16 per unit width (a² → 2 S dω).
 * The reflection coefficient R (0–1) is an argument: flexible, low-draft mats reflect little in long waves and
 * approach full reflection only for wavelengths comparable to the float size.
 * @param {object} o
 * @param {number} o.Hs significant wave height (m) (or wave height H with `regular`)
 * @param {number} o.width island width normal to the waves (m)
 * @param {number} o.reflection reflection coefficient R (-)
 * @param {boolean} [o.regular=false] treat Hs as a regular wave height H
 * @param {number} [o.rho=1000]
 * @returns {number} mean drift force (N)
 */
export function meanDriftForce({ Hs, width, reflection, regular = false, rho = 1000 }) {
  need(Hs >= 0, 'Hs must be >= 0');
  finitePos(width, 'width');
  need(reflection >= 0 && reflection <= 1, 'reflection must be within [0, 1]');
  const perWidth = regular ? 0.5 * rho * G * (Hs / 2) ** 2 : (rho * G * Hs * Hs) / 16;
  return perWidth * reflection * reflection * width;
}

/**
 * Steady current drag, F = ½ ρ C_d A U² (C_d and the reference area — submerged frontal area of floats, plus a
 * skin-friction allowance for long mats — are project inputs).
 * @param {object} o
 * @param {number} o.speed current speed (m/s)
 * @param {number} o.area reference area (m²)
 * @param {number} o.Cd drag coefficient (-)
 * @param {number} [o.rho=1000]
 * @returns {number} force (N)
 */
export function currentDrag({ speed, area, Cd, rho = 1000 }) {
  need(area >= 0 && Cd >= 0, 'area and Cd must be >= 0');
  return 0.5 * rho * Cd * area * speed * Math.abs(speed);
}

/**
 * Total horizontal environmental load on an island for one direction (wind, wave and current each with their
 * own compass direction, toward which the load acts).
 * @param {{force:number, dirDeg:number}[]} loads force (N) and direction it acts toward (deg clockwise from north)
 * @returns {{Fx:number, Fy:number, magnitude:number, dirDeg:number}} x east, y north
 */
export function combineLoads(loads) {
  need(Array.isArray(loads), 'loads must be an array');
  let Fx = 0;
  let Fy = 0;
  for (const l of loads) {
    need(Number.isFinite(l.force) && Number.isFinite(l.dirDeg), 'each load needs finite force and dirDeg');
    Fx += l.force * Math.sin(l.dirDeg * DEG);
    Fy += l.force * Math.cos(l.dirDeg * DEG);
  }
  return { Fx, Fy, magnitude: Math.hypot(Fx, Fy), dirDeg: ((Math.atan2(Fx, Fy) / DEG) % 360 + 360) % 360 };
}

/**
 * Mooring load sharing for a rigid island translating under a horizontal load. Each line i points from the
 * island toward its anchor (unit vector n_i from its compass azimuth) with horizontal stiffness k_i and
 * pretension T0_i; for a displacement d its tension is T_i = T0_i − k_i n_i·d. Equilibrium
 *   (Σ k_i n_i n_iᵀ) d = F + Σ T0_i n_i
 * is solved iteratively, lines whose tension would be negative being removed (slack). Island rotation and line
 * nonlinearity (catenary, hysteresis) are not modelled.
 * @param {object} o
 * @param {number} o.Fx load component east (N)
 * @param {number} o.Fy load component north (N)
 * @param {{azimuthDeg:number, stiffness:number, pretension?:number}[]} o.lines direction toward the anchor
 *   (deg clockwise from north), horizontal stiffness (N/m), pretension (N)
 * @returns {{offsetX:number, offsetY:number, offset:number, tensions:number[], maxTension:number,
 *   slackLines:number, loadedLines:number, sharing:number}} sharing = load magnitude / (max tension × loaded
 *   lines), 1 when all loaded lines carry the same load along the load direction
 */
export function mooringLoadSharing({ Fx, Fy, lines }) {
  need(Number.isFinite(Fx) && Number.isFinite(Fy), 'Fx and Fy must be finite');
  need(Array.isArray(lines) && lines.length >= 2, 'at least two mooring lines are required');
  const n = lines.map((l) => {
    finitePos(l.stiffness, 'line stiffness');
    return [Math.sin(l.azimuthDeg * DEG), Math.cos(l.azimuthDeg * DEG)];
  });
  const active = lines.map(() => true);
  let d = [0, 0];
  let T = [];
  for (let it = 0; it <= lines.length; it++) {
    let kxx = 0;
    let kxy = 0;
    let kyy = 0;
    let rx = Fx;
    let ry = Fy;
    for (let i = 0; i < lines.length; i++) {
      if (!active[i]) continue;
      const k = lines[i].stiffness;
      const t0 = lines[i].pretension ?? 0;
      kxx += k * n[i][0] * n[i][0];
      kxy += k * n[i][0] * n[i][1];
      kyy += k * n[i][1] * n[i][1];
      rx += t0 * n[i][0];
      ry += t0 * n[i][1];
    }
    const det = kxx * kyy - kxy * kxy;
    need(det > 1e-12 * Math.max(kxx * kyy, 1e-30), 'mooring system is a mechanism for this load (remaining lines are collinear or all slack)');
    d = [(rx * kyy - ry * kxy) / det, (ry * kxx - rx * kxy) / det];
    T = lines.map((l, i) => (active[i] ? (l.pretension ?? 0) - l.stiffness * (n[i][0] * d[0] + n[i][1] * d[1]) : 0));
    let changed = false;
    const tolT = 1e-9 * (Math.hypot(Fx, Fy) + Math.max(...lines.map((l) => l.pretension ?? 0)));
    for (let i = 0; i < lines.length; i++) {
      if (active[i] && T[i] < -tolT) {
        active[i] = false;
        changed = true;
      }
    }
    if (!changed) break;
  }
  T = T.map((t) => Math.max(0, t));
  const maxT = Math.max(...T);
  const slack = active.filter((a) => !a).length;
  const loaded = T.filter((t) => t > 1e-9 * (maxT || 1)).length;
  return {
    offsetX: d[0],
    offsetY: d[1],
    offset: Math.hypot(d[0], d[1]),
    tensions: T,
    maxTension: maxT,
    slackLines: slack,
    loadedLines: loaded,
    sharing: maxT > 0 ? Math.hypot(Fx, Fy) / (maxT * loaded) : NaN,
  };
}

/**
 * Mooring-line geometry for a varying water level (straight-line kinematics, anchor at the bed or bank at
 * horizontal distance x from the fairlead): a line that just becomes taut at the highest level has
 * L = √(x² + h_max²); at the lowest level it has geometric slack L − √(x² + h_min²), allowing a horizontal
 * excursion √(L² − h_min²) − x. An elastic line kept taut at the lowest level must instead strain by
 * ε = √(x² + h_max²)/√(x² + h_min²) − 1 at the highest level. The stiffness needed to limit the offset under a
 * load F to δ is k = F/δ.
 * @param {object} o
 * @param {number} o.anchorDistance horizontal distance fairlead–anchor x (m)
 * @param {number} o.depthMin vertical fairlead–anchor distance at the lowest water level (m)
 * @param {number} o.depthMax vertical distance at the highest water level (m)
 * @param {number} [o.load] horizontal design load per line (N)
 * @param {number} [o.allowableOffset] allowable island offset (m)
 * @returns {{length:number, slackAtLowWater:number, excursionAtLowWater:number, elasticStrain:number,
 *   requiredStiffness:number}} m, m, m, -, N/m (NaN without load and allowableOffset)
 */
export function mooringLevelRange({ anchorDistance, depthMin, depthMax, load, allowableOffset }) {
  finitePos(anchorDistance, 'anchorDistance');
  need(depthMin >= 0 && depthMax >= depthMin, 'depthMax must be >= depthMin >= 0');
  const L = Math.hypot(anchorDistance, depthMax);
  const Llow = Math.hypot(anchorDistance, depthMin);
  return {
    length: L,
    slackAtLowWater: L - Llow,
    excursionAtLowWater: Math.sqrt(L * L - depthMin * depthMin) - anchorDistance,
    elasticStrain: L / Llow - 1,
    requiredStiffness: load > 0 && allowableOffset > 0 ? load / allowableOffset : NaN,
  };
}

/**
 * Number of anchors (or mooring lines) for a horizontal design load:
 * N = ceil(F γ / (R_a η)), with anchor design capacity R_a, load factor γ and a sharing efficiency η (fraction
 * of the installed anchors effectively resisting a given load direction, from `mooringLoadSharing`).
 * @param {object} o
 * @param {number} o.load total horizontal design load (N)
 * @param {number} o.anchorCapacity design holding capacity of one anchor (N)
 * @param {number} [o.loadFactor=1] partial load factor γ (-)
 * @param {number} [o.sharing=1] sharing efficiency η (0–1]
 * @returns {{count:number, loadPerAnchor:number, utilisation:number}}
 */
export function anchorCount({ load, anchorCapacity, loadFactor = 1, sharing = 1 }) {
  need(load >= 0, 'load must be >= 0');
  finitePos(anchorCapacity, 'anchorCapacity');
  finitePos(loadFactor, 'loadFactor');
  need(sharing > 0 && sharing <= 1, 'sharing must be within (0, 1]');
  const count = Math.max(1, Math.ceil((load * loadFactor) / (anchorCapacity * sharing) - 1e-12));
  const per = (load * loadFactor) / (count * sharing);
  return { count, loadPerAnchor: per, utilisation: per / anchorCapacity };
}

/**
 * Articulated (hinged) raft of rigid modules on a regular long-crested wave, quasi-static kinematic screening.
 * Each module of length L is taken to float in its own hydrostatic equilibrium on η = a cos(kx)
 * (heave = mean of η over the module, pitch = least-squares slope); with u = kL/2:
 *   module pitch amplitude        θ = k a · 3 (sin u − u cos u)/u³,
 *   chord-slope relative rotation Δθ = (4 a/L) sin² u   (maximum 2H/L at L = λ/2),
 *   free-end mismatch at a hinge  δ = 2 a |S sin u − G cos u|,  S = sin u/u,  G = 3 (sin u − u cos u)/u².
 * A hinge transmits no moment; closing the mismatch requires a shear V = ρ g b L δ / 8 (an end force V on a
 * free prismatic module deflects its end by 4V/(ρ g b L): heave V/(ρ g b L) plus pitch 3V/(ρ g b L); two
 * modules meet at each hinge). Neighbouring-hinge interaction, dynamics, mooring and wave diffraction are
 * neglected, so the result is a screening value. Regimes by λ/L: > ~6 the raft contours the wave (small
 * connector demand); 1–6 significant articulation, largest near λ/L = 2; < 1 modules average the wave (bridging;
 * local slamming and float bending govern instead).
 * @param {object} o
 * @param {number} o.H wave height (m)
 * @param {number} o.wavelength wavelength λ (m)
 * @param {number} o.moduleLength module length along the wave direction L (m)
 * @param {number} o.moduleWidth module width b (m)
 * @param {number} [o.fillRatio=1] waterplane area / (L b) (-)
 * @param {number} [o.rho=1000]
 * @returns {{waveSlope:number, modulePitch:number, relativeRotation:number, relativeRotationDeg:number,
 *   mismatch:number, connectorShear:number, lengthRatio:number, regime:'contouring'|'articulating'|'bridging'}}
 *   rotations in rad, mismatch in m, shear in N
 */
export function hingedRaftResponse({ H, wavelength, moduleLength, moduleWidth, fillRatio = 1, rho = 1000 }) {
  need(H >= 0, 'H must be >= 0');
  finitePos(wavelength, 'wavelength');
  finitePos(moduleLength, 'moduleLength');
  finitePos(moduleWidth, 'moduleWidth');
  need(fillRatio > 0 && fillRatio <= 1, 'fillRatio must be within (0, 1]');
  const a = H / 2;
  const k = (2 * Math.PI) / wavelength;
  const u = (k * moduleLength) / 2;
  const S = Math.sin(u) / u;
  const Gf = (3 * (Math.sin(u) - u * Math.cos(u))) / (u * u);
  const rel = ((4 * a) / moduleLength) * Math.sin(u) ** 2;
  const mismatch = 2 * a * Math.abs(S * Math.sin(u) - Gf * Math.cos(u));
  const ratio = wavelength / moduleLength;
  return {
    waveSlope: k * a,
    modulePitch: (k * a * Gf) / u,
    relativeRotation: rel,
    relativeRotationDeg: rel / DEG,
    mismatch,
    connectorShear: (rho * G * moduleWidth * fillRatio * moduleLength * mismatch) / 8,
    lengthRatio: ratio,
    regime: ratio > 6 ? 'contouring' : ratio >= 1 ? 'articulating' : 'bridging',
  };
}

/**
 * Water saving from reduced evaporation under a floating array:
 * V = E A_covered f, with the open-water evaporation rate E supplied by the user (e.g. Penman 1948 estimate or
 * pan data) and f the fractional reduction under the covered area (float-type dependent; an argument).
 * @param {object} o
 * @param {number} o.evaporationRate open-water evaporation (m/yr)
 * @param {number} o.coveredArea water surface covered by floats (m²)
 * @param {number} o.reduction fractional evaporation reduction under the cover (0–1)
 * @returns {number} water saved (m³/yr)
 */
export function evaporationSaving({ evaporationRate, coveredArea, reduction }) {
  need(evaporationRate >= 0 && coveredArea >= 0, 'evaporationRate and coveredArea must be >= 0');
  need(reduction >= 0 && reduction <= 1, 'reduction must be within [0, 1]');
  return evaporationRate * coveredArea * reduction;
}

/**
 * Yield gain from lower operating temperature over water: ΔP/P = γ ΔT_c with ΔT_c the (negative) difference in
 * irradiance-weighted cell temperature between the floating and the reference land system (from the thermal
 * models in pv.js with the respective coefficients).
 * @param {object} o
 * @param {number} o.deltaTempC floating minus land cell temperature (°C), normally negative
 * @param {number} [o.gammaPdc=-0.0035] power temperature coefficient (1/K)
 * @returns {number} fractional energy gain (-)
 */
export function coolingYieldGain({ deltaTempC, gammaPdc = -0.0035 }) {
  need(Number.isFinite(deltaTempC), 'deltaTempC must be finite');
  return gammaPdc * deltaTempC;
}

/**
 * Length of a floating/dynamic cable between the island and a fixed point, allowing for the island excursion
 * and the water-level range: L = slackFactor × √((x + offset)² + ΔH²).
 * @param {object} o
 * @param {number} o.horizontalDistance nominal horizontal distance (m)
 * @param {number} [o.offset=0] island excursion away from the fixed point (m)
 * @param {number} [o.levelRange=0] vertical travel of the island relative to the fixed point (m)
 * @param {number} [o.slackFactor=1.1] installation slack allowance (-), >= 1
 * @returns {{length:number, slack:number}} cable length (m) and the excess over the nominal distance (m)
 */
export function floatingCableLength({ horizontalDistance, offset = 0, levelRange = 0, slackFactor = 1.1 }) {
  finitePos(horizontalDistance, 'horizontalDistance');
  need(offset >= 0 && levelRange >= 0, 'offset and levelRange must be >= 0');
  need(slackFactor >= 1, 'slackFactor must be >= 1');
  const length = slackFactor * Math.hypot(horizontalDistance + offset, levelRange);
  return { length, slack: length - horizontalDistance };
}

/**
 * Site-exposure screening by significant wave height. Thresholds are arguments; the defaults (1 m and 3 m)
 * follow the common grouping into sheltered inland waters, nearshore/coastal and offshore sites.
 * @param {object} o
 * @param {number} o.Hs design significant wave height (m)
 * @param {number[]} [o.thresholds=[1, 3]] ascending class boundaries (m)
 * @param {string[]} [o.labels=['inland', 'nearshore', 'offshore']] one more label than thresholds
 * @returns {{class:string, index:number}}
 */
export function exposureClass({ Hs, thresholds = [1, 3], labels = ['inland', 'nearshore', 'offshore'] }) {
  need(Hs >= 0, 'Hs must be >= 0');
  need(labels.length === thresholds.length + 1, 'labels must have one more entry than thresholds');
  let i = 0;
  while (i < thresholds.length && Hs >= thresholds[i]) i++;
  return { class: labels[i], index: i };
}

/**
 * Electrical design reminders specific to floating arrays (qualitative; no numerical criteria).
 */
export const FLOATING_PV_ELECTRICAL_NOTES = Object.freeze([
  'Size floating and dynamic cables for the full island excursion and water-level range plus installation slack; avoid tension at connectors.',
  'Use equipotential bonding across articulated float connections with flexible bonding jumpers; hinge hardware is not a reliable earth path.',
  'Check insulation resistance and leakage current for modules and cables operating in a permanently humid, splash-exposed environment.',
  'Route DC strings so that cables are not submerged unless they are rated for continuous immersion.',
]);
