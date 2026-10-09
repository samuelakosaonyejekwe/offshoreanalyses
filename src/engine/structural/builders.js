/**
 * Parametric generators of frame3d models for common offshore and marine structures: 4-leg and 3-leg jackets,
 * monopile + tower stick models, deck / topsides grid frames, wharf or jetty bents, dolphins, portal frames, and
 * the conversion of an imported `frame` geometry model.
 *
 * Every generator returns a model in the frame3d format ({ nodes, members, supports, nodeMasses?, loadCases })
 * plus a `meta` object describing the generated topology (node indices by level, member groups). Z is up; for
 * the marine structures z = 0 is the still-water level and the seabed is at z = −waterDepth. SI units.
 * The generated models are starting points for screening analyses: joint cans, conductor framing, appurtenances
 * and local eccentricities are not generated.
 */

const STEEL = Object.freeze({ E: 2.1e11, nu: 0.3, rho: 7850, alpha: 1.2e-5, fy: 355e6 });

function pos(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}

function tubeSec(name, s) {
  if (!s || !(s.D > 0) || !(s.t > 0)) throw new RangeError(`${name} must be { D, t } with positive values`);
  return { type: 'tube', D: s.D, t: s.t };
}

class Builder {
  constructor() { this.nodes = []; this.members = []; this.index = new Map(); }

  node(x, y, z) {
    const key = `${Math.round(x * 1e5)},${Math.round(y * 1e5)},${Math.round(z * 1e5)}`;
    let i = this.index.get(key);
    if (i === undefined) { i = this.nodes.length; this.nodes.push([x, y, z]); this.index.set(key, i); }
    return i;
  }

  member(a, b, section, group, extra) {
    if (a === b) return -1;
    this.members.push({ a, b, section, group, ...extra });
    return this.members.length - 1;
  }
}

function pileSupport(b, supports, node, piles, material) {
  const mode = (piles && piles.mode) || 'fixed';
  if (mode === 'fixed') { supports.push({ node, fix: 'fixed' }); return; }
  if (mode === 'pinned') { supports.push({ node, fix: 'pinned' }); return; }
  if (mode === 'springs') {
    if (!piles.spring) throw new RangeError('piles.spring (6 diagonal values or a 6×6 matrix) is required for mode "springs"');
    supports.push({ node, spring: piles.spring });
    return;
  }
  if (mode === 'fixity') {
    const sec = tubeSec('piles.section', piles.section);
    const depth = piles.fixityDepth !== undefined ? pos('piles.fixityDepth', piles.fixityDepth) : 6 * sec.D;
    const p = b.nodes[node];
    const tip = b.node(p[0], p[1], p[2] - depth);
    b.member(tip, node, sec, 'pile', { material, flooded: true });
    supports.push({ node: tip, fix: 'fixed' });
    return;
  }
  throw new RangeError(`piles.mode "${mode}" is not one of fixed, pinned, springs, fixity`);
}

/**
 * Jacket substructure with 4 legs (rectangular plan) or 3 legs (equilateral plan), battered legs, horizontal
 * framing at every bay level and X, K or single-diagonal face bracing.
 * Pile modelling: 'fixed' / 'pinned' supports at the mudline, 'springs' (pile-head stiffness, 6 diagonal values or a
 * full 6×6 matrix per leg, as for piles grouted in the legs and linearised), or 'fixity' (explicit pile stubs down
 * to an apparent fixity depth, default 6 pile diameters, fixed at the tip).
 * @param {object} o
 * @param {3|4} [o.legs=4] number of legs
 * @param {number} o.waterDepth still-water depth [m]
 * @param {number} o.topElevation elevation of the jacket top (work point) above still water [m]
 * @param {number} o.topWidthX leg spacing at the top in x [m] (for 3 legs: side length of the triangle at the top)
 * @param {number} [o.topWidthY=o.topWidthX] leg spacing at the top in y [m] (4 legs)
 * @param {number} [o.batter=8] leg batter 1 horizontal : batter vertical in each battered direction [-] (Infinity = vertical legs)
 * @param {number} [o.baseWidthX] leg spacing at the mudline in x [m]; overrides `batter` when given
 * @param {number} [o.baseWidthY] leg spacing at the mudline in y [m]
 * @param {number} [o.bays=4] number of bracing bays
 * @param {'X'|'K'|'diagonal'} [o.bracing='X'] face bracing pattern
 * @param {boolean} [o.planBracing=true] one plan diagonal per horizontal level (4 legs)
 * @param {{D:number, t:number}} o.legSection [m]
 * @param {{D:number, t:number}} o.braceSection [m]
 * @param {{D:number, t:number}} [o.horizontalSection=o.braceSection] [m]
 * @param {{mode?:string, spring?:any, section?:object, fixityDepth?:number}} [o.piles] pile / foundation idealisation
 * @param {object} [o.material] steel properties (default S355: E 210 GPa, ρ 7850 kg/m³)
 * @param {number} [o.deckMass=0] topsides mass shared equally by the leg tops [kg]
 * @param {boolean} [o.floodedLegs=true] legs flooded below the waterline
 * @returns {{nodes:number[][], members:object[], supports:object[], nodeMasses:object[], loadCases:object[], meta:object}} frame3d model;
 *   meta = { levels (z of each framing level, bottom first), legNodes[level][leg], topNodes, mudlineNodes, groups }
 */
export function jacketModel({
  legs = 4, waterDepth, topElevation, topWidthX, topWidthY, batter = 8, baseWidthX, baseWidthY, bays = 4, bracing = 'X',
  planBracing = true, legSection, braceSection, horizontalSection, piles, material = STEEL, deckMass = 0, floodedLegs = true,
}) {
  if (legs !== 3 && legs !== 4) throw new RangeError('legs must be 3 or 4');
  pos('waterDepth', waterDepth); pos('topElevation', topElevation); pos('topWidthX', topWidthX);
  if (!Number.isInteger(bays) || bays < 1) throw new RangeError('bays must be a positive integer');
  if (!['X', 'K', 'diagonal'].includes(bracing)) throw new RangeError('bracing must be "X", "K" or "diagonal"');
  if (!(batter > 0)) throw new RangeError('batter must be positive');
  const leg = tubeSec('legSection', legSection);
  const brace = tubeSec('braceSection', braceSection);
  const horiz = horizontalSection ? tubeSec('horizontalSection', horizontalSection) : brace;
  const H = waterDepth + topElevation;
  const tw = [topWidthX, topWidthY === undefined ? topWidthX : pos('topWidthY', topWidthY)];
  const grow = Number.isFinite(batter) ? 2 * H / batter : 0;
  const bw = [
    baseWidthX === undefined ? tw[0] + grow : pos('baseWidthX', baseWidthX),
    baseWidthY === undefined ? tw[1] + grow : pos('baseWidthY', baseWidthY),
  ];
  const b = new Builder();
  const levels = [];
  const legNodes = [];
  // plan position of leg k at relative height s (0 = mudline, 1 = top)
  const plan = (k, s) => {
    const wx = bw[0] + s * (tw[0] - bw[0]);
    if (legs === 4) {
      const wy = bw[1] + s * (tw[1] - bw[1]);
      const sx = [-1, 1, 1, -1][k]; const sy = [-1, -1, 1, 1][k];
      return [sx * wx / 2, sy * wy / 2];
    }
    // equilateral triangle of side wx centred on the origin
    const R = wx / Math.sqrt(3);
    const a = (90 + 120 * k) * Math.PI / 180;
    return [R * Math.cos(a), R * Math.sin(a)];
  };
  for (let l = 0; l <= bays; l++) {
    const s = l / bays;
    const z = -waterDepth + s * H;
    levels.push(z);
    const row = [];
    for (let k = 0; k < legs; k++) { const [x, y] = plan(k, s); row.push(b.node(x, y, z)); }
    legNodes.push(row);
  }
  const mat = { material };
  for (let l = 0; l < bays; l++) {
    for (let k = 0; k < legs; k++) b.member(legNodes[l][k], legNodes[l + 1][k], leg, 'leg', { ...mat, flooded: floodedLegs });
  }
  // K bracing needs a node at the middle of the upper horizontal of each bay
  const midNode = (l, k) => {
    const p = b.nodes[legNodes[l][k]]; const q = b.nodes[legNodes[l][(k + 1) % legs]];
    return b.node((p[0] + q[0]) / 2, (p[1] + q[1]) / 2, (p[2] + q[2]) / 2);
  };
  for (let l = 0; l <= bays; l++) {
    for (let k = 0; k < legs; k++) {
      const n1 = legNodes[l][k]; const n2 = legNodes[l][(k + 1) % legs];
      if (bracing === 'K' && l > 0) {
        const m = midNode(l, k);
        b.member(n1, m, horiz, 'horizontal', mat); b.member(m, n2, horiz, 'horizontal', mat);
      } else b.member(n1, n2, horiz, 'horizontal', mat);
    }
    if (planBracing && legs === 4) b.member(legNodes[l][l % 2], legNodes[l][l % 2 + 2], horiz, 'plan', mat);
  }
  for (let l = 0; l < bays; l++) {
    for (let k = 0; k < legs; k++) {
      const a1 = legNodes[l][k]; const a2 = legNodes[l][(k + 1) % legs];
      const c1 = legNodes[l + 1][k]; const c2 = legNodes[l + 1][(k + 1) % legs];
      if (bracing === 'X') { b.member(a1, c2, brace, 'brace', mat); b.member(a2, c1, brace, 'brace', mat); }
      else if (bracing === 'K') { const m = midNode(l + 1, k); b.member(a1, m, brace, 'brace', mat); b.member(a2, m, brace, 'brace', mat); }
      else if ((l + k) % 2 === 0) b.member(a1, c2, brace, 'brace', mat);
      else b.member(a2, c1, brace, 'brace', mat);
    }
  }
  const supports = [];
  const mudlineNodes = legNodes[0].slice();
  for (const n of mudlineNodes) pileSupport(b, supports, n, piles, material);
  const nodeMasses = [];
  if (deckMass > 0) for (const n of legNodes[bays]) nodeMasses.push({ node: n, m: deckMass / legs });
  const groups = {};
  b.members.forEach((m, i) => { (groups[m.group] = groups[m.group] || []).push(i); });
  return {
    nodes: b.nodes, members: b.members, supports, nodeMasses,
    loadCases: [{ name: 'selfWeight', selfWeight: { buoyancy: { waterline: 0 } } }],
    meta: { kind: 'jacket', legs, levels, legNodes, topNodes: legNodes[bays].slice(), mudlineNodes, groups, baseWidth: bw, topWidth: tw },
  };
}

/**
 * Monopile + tower stick model for an offshore wind turbine or a single-pile structure: prismatic pile from the
 * foundation level to the interface, linearly tapered tower represented by prismatic segments with the mean
 * diameter and thickness of each segment, lumped rotor–nacelle mass at the top.
 * Foundation: 'fixity' (pile extended below the mudline to an apparent fixity depth and clamped), 'spring'
 * (6×6 or diagonal pile-head stiffness at the mudline) or 'fixed' at the mudline.
 * @param {object} o
 * @param {number} o.waterDepth [m]
 * @param {number} o.pileD monopile outer diameter [m]
 * @param {number} o.pileT monopile wall thickness [m]
 * @param {number} o.interfaceZ elevation of the tower base above still water [m]
 * @param {number} o.hubHeight elevation of the tower top / RNA mass above still water [m]
 * @param {number} o.towerBaseD tower base diameter [m]
 * @param {number} o.towerTopD tower top diameter [m]
 * @param {number} o.towerBaseT tower base wall thickness [m]
 * @param {number} [o.towerTopT=o.towerBaseT] tower top wall thickness [m]
 * @param {number} [o.rnaMass=0] rotor–nacelle assembly mass [kg]
 * @param {number[]} [o.rnaInertia] RNA rotary inertia [Ixx, Iyy, Izz] [kg·m²]
 * @param {{mode?:'fixity'|'spring'|'fixed', fixityDepth?:number, spring?:any}} [o.foundation] default fixity at 4 pile diameters
 * @param {number} [o.nPile=8] elements in the pile between the mudline and the interface
 * @param {number} [o.nTower=12] elements in the tower
 * @param {object} [o.material]
 * @returns {{nodes:number[][], members:object[], supports:object[], nodeMasses:object[], loadCases:object[], meta:object}} frame3d model;
 *   meta = { mudlineNode, interfaceNode, topNode, pileMembers, towerMembers }
 */
export function monopileTowerModel({
  waterDepth, pileD, pileT, interfaceZ, hubHeight, towerBaseD, towerTopD, towerBaseT, towerTopT, rnaMass = 0, rnaInertia,
  foundation = {}, nPile = 8, nTower = 12, material = STEEL,
}) {
  pos('waterDepth', waterDepth); pos('pileD', pileD); pos('pileT', pileT); pos('hubHeight', hubHeight);
  pos('towerBaseD', towerBaseD); pos('towerTopD', towerTopD); pos('towerBaseT', towerBaseT);
  if (!(interfaceZ > -waterDepth && interfaceZ < hubHeight)) throw new RangeError('interfaceZ must lie between the mudline and hubHeight');
  if (!Number.isInteger(nPile) || nPile < 1 || !Number.isInteger(nTower) || nTower < 1) throw new RangeError('nPile and nTower must be positive integers');
  const tTop = towerTopT === undefined ? towerBaseT : pos('towerTopT', towerTopT);
  const b = new Builder();
  const supports = [];
  const pileSec = { type: 'tube', D: pileD, t: pileT };
  const mode = foundation.mode || 'fixity';
  const mud = b.node(0, 0, -waterDepth);
  const pileMembers = []; const towerMembers = [];
  if (mode === 'fixity') {
    const depth = foundation.fixityDepth === undefined ? 4 * pileD : pos('foundation.fixityDepth', foundation.fixityDepth);
    const tip = b.node(0, 0, -waterDepth - depth);
    pileMembers.push(b.member(tip, mud, pileSec, 'pile', { material, flooded: true }));
    supports.push({ node: tip, fix: 'fixed' });
  } else if (mode === 'spring') {
    if (!foundation.spring) throw new RangeError('foundation.spring is required for mode "spring"');
    supports.push({ node: mud, spring: foundation.spring });
  } else if (mode === 'fixed') supports.push({ node: mud, fix: 'fixed' });
  else throw new RangeError(`foundation.mode "${mode}" is not one of fixity, spring, fixed`);
  let prev = mud;
  for (let i = 1; i <= nPile; i++) {
    const n = b.node(0, 0, -waterDepth + (interfaceZ + waterDepth) * i / nPile);
    pileMembers.push(b.member(prev, n, pileSec, 'pile', { material, flooded: true }));
    prev = n;
  }
  const iface = prev;
  for (let i = 1; i <= nTower; i++) {
    const s0 = (i - 1) / nTower; const s1 = i / nTower; const sm = 0.5 * (s0 + s1);
    const n = b.node(0, 0, interfaceZ + (hubHeight - interfaceZ) * s1);
    const sec = { type: 'tube', D: towerBaseD + sm * (towerTopD - towerBaseD), t: towerBaseT + sm * (tTop - towerBaseT) };
    towerMembers.push(b.member(prev, n, sec, 'tower', { material }));
    prev = n;
  }
  const nodeMasses = [];
  if (rnaMass > 0) nodeMasses.push({ node: prev, m: rnaMass, inertia: rnaInertia });
  return {
    nodes: b.nodes, members: b.members, supports, nodeMasses,
    loadCases: [{ name: 'selfWeight', selfWeight: { buoyancy: { waterline: 0 } } }],
    meta: { kind: 'monopile', mudlineNode: mud, interfaceNode: iface, topNode: prev, pileMembers, towerMembers },
  };
}

/**
 * Deck / topsides grid frame: an orthogonal grid of beams at one elevation (girders along x, beams along y) on
 * columns down to the support level.
 * @param {object} o
 * @param {number} o.lengthX plan length [m]
 * @param {number} o.widthY plan width [m]
 * @param {number} [o.nx=3] bays along x
 * @param {number} [o.ny=2] bays along y
 * @param {number} [o.elevation=0] deck elevation z [m]
 * @param {number} o.columnHeight column height below the deck [m]
 * @param {object} o.girderSection section of the members along x (frame3d section description)
 * @param {object} [o.beamSection=o.girderSection] section of the members along y
 * @param {object} o.columnSection column section
 * @param {'corners'|'perimeter'|'all'} [o.columns='corners'] grid points that receive a column
 * @param {'fixed'|'pinned'} [o.columnBase='fixed']
 * @param {number} [o.areaLoad=0] uniform deck pressure for the 'deckLoad' case [Pa], distributed to the grid nodes by tributary area
 * @param {object} [o.material]
 * @returns {{nodes:number[][], members:object[], supports:object[], loadCases:object[], meta:object}} frame3d model; meta.grid[i][j] = node index
 */
export function deckGridModel({
  lengthX, widthY, nx = 3, ny = 2, elevation = 0, columnHeight, girderSection, beamSection, columnSection, columns = 'corners',
  columnBase = 'fixed', areaLoad = 0, material = STEEL,
}) {
  pos('lengthX', lengthX); pos('widthY', widthY); pos('columnHeight', columnHeight);
  if (!Number.isInteger(nx) || nx < 1 || !Number.isInteger(ny) || ny < 1) throw new RangeError('nx and ny must be positive integers');
  if (!girderSection || !columnSection) throw new RangeError('girderSection and columnSection are required');
  if (!['corners', 'perimeter', 'all'].includes(columns)) throw new RangeError('columns must be "corners", "perimeter" or "all"');
  if (columnBase !== 'fixed' && columnBase !== 'pinned') throw new RangeError('columnBase must be "fixed" or "pinned"');
  const b = new Builder();
  const grid = [];
  const dx = lengthX / nx; const dy = widthY / ny;
  for (let i = 0; i <= nx; i++) { grid.push([]); for (let j = 0; j <= ny; j++) grid[i].push(b.node(i * dx, j * dy, elevation)); }
  for (let j = 0; j <= ny; j++) for (let i = 0; i < nx; i++) b.member(grid[i][j], grid[i + 1][j], girderSection, 'girder', { material });
  for (let i = 0; i <= nx; i++) for (let j = 0; j < ny; j++) b.member(grid[i][j], grid[i][j + 1], beamSection || girderSection, 'beam', { material });
  const supports = []; const columnNodes = [];
  for (let i = 0; i <= nx; i++) {
    for (let j = 0; j <= ny; j++) {
      const edgeI = i === 0 || i === nx; const edgeJ = j === 0 || j === ny;
      const take = columns === 'all' || (columns === 'perimeter' && (edgeI || edgeJ)) || (columns === 'corners' && edgeI && edgeJ);
      if (!take) continue;
      const foot = b.node(i * dx, j * dy, elevation - columnHeight);
      b.member(foot, grid[i][j], columnSection, 'column', { material });
      supports.push({ node: foot, fix: columnBase });
      columnNodes.push(grid[i][j]);
    }
  }
  const loadCases = [{ name: 'selfWeight', selfWeight: {} }];
  if (areaLoad !== 0) {
    const nodal = [];
    for (let i = 0; i <= nx; i++) {
      for (let j = 0; j <= ny; j++) {
        const ax = (i === 0 || i === nx ? 0.5 : 1) * dx; const ay = (j === 0 || j === ny ? 0.5 : 1) * dy;
        nodal.push({ node: grid[i][j], F: [0, 0, -areaLoad * ax * ay, 0, 0, 0] });
      }
    }
    loadCases.push({ name: 'deckLoad', nodal });
  }
  return { nodes: b.nodes, members: b.members, supports, loadCases, meta: { kind: 'deck', grid, columnNodes } };
}

/**
 * Wharf / jetty bent in the x–z plane: a deck beam on vertical piles with optional raker (batter) piles. Piles are
 * fixed at an apparent fixity depth below the seabed at their location; out-of-plane displacement and twist are
 * restrained at the deck so that the bent can be analysed on its own.
 * @param {object} o
 * @param {number[]} o.pileX x positions of the pile heads along the deck beam [m] (ascending)
 * @param {number} o.deckElevation z of the deck beam axis [m]
 * @param {number|number[]} o.seabedZ seabed elevation at each pile (or one value) [m]
 * @param {number} [o.fixityDepth] depth of apparent fixity below the seabed [m] (default 5 pile diameters)
 * @param {Array<{pile:number, batter:number, direction:(1|-1)}>} [o.rakers] raker piles: index into pileX, batter 1 horizontal : batter
 *   vertical, direction of the toe in x
 * @param {object} o.deckSection deck beam section
 * @param {{D:number, t:number}} o.pileSection pile section [m]
 * @param {number} [o.overhang=0] deck overhang beyond the end piles [m]
 * @param {boolean} [o.pinnedHeads=false] pile heads pinned to the deck beam
 * @param {object} [o.material]
 * @returns {{nodes:number[][], members:object[], supports:object[], loadCases:object[], meta:object}} frame3d model; meta = { headNodes, toeNodes, deckMembers, pileMembers }
 */
export function wharfBentModel({
  pileX, deckElevation, seabedZ, fixityDepth, rakers = [], deckSection, pileSection, overhang = 0, pinnedHeads = false, material = STEEL,
}) {
  if (!Array.isArray(pileX) || pileX.length < 1) throw new RangeError('pileX must list at least one pile position');
  for (let i = 1; i < pileX.length; i++) if (!(pileX[i] > pileX[i - 1])) throw new RangeError('pileX must be strictly ascending');
  if (!Number.isFinite(deckElevation)) throw new RangeError('deckElevation must be finite');
  if (!deckSection) throw new RangeError('deckSection is required');
  const pile = tubeSec('pileSection', pileSection);
  const fix = fixityDepth === undefined ? 5 * pile.D : pos('fixityDepth', fixityDepth);
  const bedAt = (i) => (Array.isArray(seabedZ) ? seabedZ[i] : seabedZ);
  const b = new Builder();
  const headNodes = pileX.map((x) => b.node(x, 0, deckElevation));
  const deckMembers = []; const pileMembers = []; const toeNodes = []; const supports = [];
  const deckNodes = headNodes.slice();
  if (overhang > 0) { deckNodes.unshift(b.node(pileX[0] - overhang, 0, deckElevation)); deckNodes.push(b.node(pileX.at(-1) + overhang, 0, deckElevation)); }
  for (let i = 0; i + 1 < deckNodes.length; i++) deckMembers.push(b.member(deckNodes[i], deckNodes[i + 1], deckSection, 'deck', { material }));
  const rel = pinnedHeads ? { releases: { b: 'ry' } } : {};
  const addPile = (head, x, bed, slope, group) => {
    if (!Number.isFinite(bed) || !(bed < deckElevation)) throw new RangeError('seabedZ must be below the deck at every pile');
    const zt = bed - fix;
    const toe = b.node(x + slope * (deckElevation - zt), 0, zt);
    pileMembers.push(b.member(toe, head, pile, group, { material, ...rel }));
    supports.push({ node: toe, fix: 'fixed' });
    toeNodes.push(toe);
  };
  pileX.forEach((x, i) => addPile(headNodes[i], x, bedAt(i), 0, 'pile'));
  for (const [k, r] of rakers.entries()) {
    if (!Number.isInteger(r.pile) || r.pile < 0 || r.pile >= pileX.length) throw new RangeError(`rakers[${k}].pile is out of range`);
    if (!(r.batter > 0)) throw new RangeError(`rakers[${k}].batter must be positive`);
    addPile(headNodes[r.pile], pileX[r.pile], bedAt(r.pile), (r.direction === -1 ? -1 : 1) / r.batter, 'raker');
  }
  for (const n of deckNodes) supports.push({ node: n, fix: [0, 1, 0, 1, 0, 1] });
  return {
    nodes: b.nodes, members: b.members, supports, loadCases: [{ name: 'selfWeight', selfWeight: { buoyancy: { waterline: 0 } } }],
    meta: { kind: 'wharfBent', headNodes, toeNodes, deckNodes, deckMembers, pileMembers },
  };
}

/**
 * Berthing / mooring dolphin: a ring of piles (optionally raked outward) under a rigid cap. The cap is represented
 * by stiff members from a central reference node to each pile head; loads are applied at the central node.
 * @param {object} o
 * @param {number} o.nPiles number of piles (>= 3)
 * @param {number} o.radius radius of the pile-head circle [m]
 * @param {number} o.capElevation z of the pile heads / cap reference node [m]
 * @param {number} o.seabedZ seabed elevation [m]
 * @param {number} [o.fixityDepth] apparent fixity depth below the seabed [m] (default 5 pile diameters)
 * @param {number} [o.batter=Infinity] pile rake 1 horizontal : batter vertical, toes outward
 * @param {{D:number, t:number}} o.pileSection [m]
 * @param {number} [o.capStiffnessFactor=1000] multiplier on the pile section properties for the cap members [-]
 * @param {number} [o.capMass=0] mass of the cap at the central node [kg]
 * @param {object} [o.material]
 * @returns {{nodes:number[][], members:object[], supports:object[], nodeMasses:object[], loadCases:object[], meta:object}} frame3d model; meta = { capNode, headNodes, toeNodes }
 */
export function dolphinModel({
  nPiles, radius, capElevation, seabedZ, fixityDepth, batter = Infinity, pileSection, capStiffnessFactor = 1000, capMass = 0, material = STEEL,
}) {
  if (!Number.isInteger(nPiles) || nPiles < 3) throw new RangeError('nPiles must be an integer of at least 3');
  pos('radius', radius);
  if (!(seabedZ < capElevation)) throw new RangeError('seabedZ must be below capElevation');
  if (!(batter > 0)) throw new RangeError('batter must be positive');
  const pile = tubeSec('pileSection', pileSection);
  const fix = fixityDepth === undefined ? 5 * pile.D : pos('fixityDepth', fixityDepth);
  const b = new Builder();
  const capNode = b.node(0, 0, capElevation);
  const di = pile.D - 2 * pile.t;
  const A = Math.PI / 4 * (pile.D ** 2 - di ** 2); const I = Math.PI / 64 * (pile.D ** 4 - di ** 4);
  const k = capStiffnessFactor;
  const capSec = { A: A * k, Iy: I * k, Iz: I * k, J: 2 * I * k };
  const massless = { ...material, rho: 0 };
  const headNodes = []; const toeNodes = []; const supports = [];
  const zt = seabedZ - fix;
  const slope = Number.isFinite(batter) ? 1 / batter : 0;
  for (let i = 0; i < nPiles; i++) {
    const a = 2 * Math.PI * i / nPiles;
    const head = b.node(radius * Math.cos(a), radius * Math.sin(a), capElevation);
    const rt = radius + slope * (capElevation - zt);
    const toe = b.node(rt * Math.cos(a), rt * Math.sin(a), zt);
    b.member(capNode, head, capSec, 'cap', { material: massless });
    b.member(toe, head, pile, 'pile', { material });
    supports.push({ node: toe, fix: 'fixed' });
    headNodes.push(head); toeNodes.push(toe);
  }
  for (let i = 0; i < nPiles; i++) b.member(headNodes[i], headNodes[(i + 1) % nPiles], capSec, 'cap', { material: massless });
  return {
    nodes: b.nodes, members: b.members, supports, nodeMasses: capMass > 0 ? [{ node: capNode, m: capMass }] : [],
    loadCases: [{ name: 'selfWeight', selfWeight: { buoyancy: { waterline: 0 } } }],
    meta: { kind: 'dolphin', capNode, headNodes, toeNodes },
  };
}

/**
 * Single- or multi-bay portal frame in the x–z plane with optional duo-pitch rafters. Out-of-plane translation
 * and the out-of-plane rotations are restrained at every node (plane-frame behaviour).
 * @param {object} o
 * @param {number} o.span bay width [m]
 * @param {number} o.height eaves height [m]
 * @param {number} [o.rise=0] apex rise above the eaves [m] (0 → flat beam)
 * @param {number} [o.bays=1] number of bays
 * @param {object} o.columnSection column section
 * @param {object} o.rafterSection rafter / beam section
 * @param {'fixed'|'pinned'} [o.base='fixed']
 * @param {number} [o.nRafter=2] elements per rafter slope (>= 1)
 * @param {object} [o.material]
 * @returns {{nodes:number[][], members:object[], supports:object[], loadCases:object[], meta:object}} frame3d model; meta = { baseNodes, eavesNodes, apexNodes, columnMembers, rafterMembers }
 */
export function portalFrameModel({ span, height, rise = 0, bays = 1, columnSection, rafterSection, base = 'fixed', nRafter = 2, material = STEEL }) {
  pos('span', span); pos('height', height);
  if (!(rise >= 0)) throw new RangeError('rise must be non-negative');
  if (!Number.isInteger(bays) || bays < 1) throw new RangeError('bays must be a positive integer');
  if (!Number.isInteger(nRafter) || nRafter < 1) throw new RangeError('nRafter must be a positive integer');
  if (base !== 'fixed' && base !== 'pinned') throw new RangeError('base must be "fixed" or "pinned"');
  if (!columnSection || !rafterSection) throw new RangeError('columnSection and rafterSection are required');
  const b = new Builder();
  const baseNodes = []; const eavesNodes = []; const apexNodes = []; const columnMembers = []; const rafterMembers = [];
  const supports = [];
  for (let c = 0; c <= bays; c++) {
    const foot = b.node(c * span, 0, 0); const top = b.node(c * span, 0, height);
    baseNodes.push(foot); eavesNodes.push(top);
    columnMembers.push(b.member(foot, top, columnSection, 'column', { material }));
    supports.push({ node: foot, fix: base === 'fixed' ? 'fixed' : [1, 1, 1, 1, 0, 1] });
  }
  for (let c = 0; c < bays; c++) {
    const x0 = c * span;
    const n = 2 * nRafter;
    let prev = eavesNodes[c];
    for (let i = 1; i <= n; i++) {
      const s = i / n;
      const z = height + rise * (1 - Math.abs(2 * s - 1));
      const node = i === n ? eavesNodes[c + 1] : b.node(x0 + s * span, 0, z);
      rafterMembers.push(b.member(prev, node, rafterSection, 'rafter', { material }));
      if (i === nRafter) apexNodes.push(node);
      prev = node;
    }
  }
  const held = new Set(baseNodes);
  for (let i = 0; i < b.nodes.length; i++) if (!held.has(i)) supports.push({ node: i, fix: [0, 1, 0, 1, 0, 1] });
  return {
    nodes: b.nodes, members: b.members, supports, loadCases: [{ name: 'selfWeight', selfWeight: {} }],
    meta: { kind: 'portal', baseNodes, eavesNodes, apexNodes, columnMembers, rafterMembers },
  };
}

/**
 * Convert an imported geometry model of kind 'frame' ({ kind, name, nodes: Float64Array(3n), members: [{ a, b,
 * section?, material?, group? }] }) into a frame3d model. Sections given as { type: 'tube', D, t } or as
 * property records { A, Iy, Iz, J } are passed through; members without a section receive `defaultSection`.
 * Supports: by default every node within `supportTolerance` of the lowest elevation is fixed.
 * @param {object} frame imported model of kind 'frame'
 * @param {object} [o]
 * @param {object} [o.defaultSection={type:'tube', D:0.5, t:0.02}] section for members without one
 * @param {object} [o.defaultMaterial] material for members without one
 * @param {Array<{node:number, fix?:any, spring?:any}>|'lowest'} [o.supports='lowest'] explicit supports, or automatic at the lowest level
 * @param {number} [o.supportTolerance] elevation tolerance for the automatic supports [m] (default 1e-3 of the model height, at least 1 mm)
 * @param {'fixed'|'pinned'} [o.supportType='fixed']
 * @returns {{nodes:number[][], members:object[], supports:object[], loadCases:object[], meta:object}} frame3d model
 */
export function fromImportedFrame(frame, { defaultSection = { type: 'tube', D: 0.5, t: 0.02 }, defaultMaterial, supports = 'lowest', supportTolerance, supportType = 'fixed' } = {}) {
  if (!frame || frame.kind !== 'frame') throw new RangeError('frame must be an imported model of kind "frame"');
  const raw = frame.nodes;
  if (!raw || raw.length < 6 || raw.length % 3 !== 0) throw new RangeError('frame.nodes must hold 3 coordinates per node for at least two nodes');
  if (!Array.isArray(frame.members) || frame.members.length === 0) throw new RangeError('frame.members must not be empty');
  const n = raw.length / 3;
  const nodes = [];
  let zmin = Infinity; let zmax = -Infinity;
  for (let i = 0; i < n; i++) {
    nodes.push([raw[3 * i], raw[3 * i + 1], raw[3 * i + 2]]);
    zmin = Math.min(zmin, raw[3 * i + 2]); zmax = Math.max(zmax, raw[3 * i + 2]);
  }
  const used = new Uint8Array(n);
  const members = [];
  frame.members.forEach((m, i) => {
    if (!Number.isInteger(m.a) || !Number.isInteger(m.b) || m.a < 0 || m.b < 0 || m.a >= n || m.b >= n) throw new RangeError(`frame.members[${i}] has invalid node indices`);
    if (m.a === m.b) return;
    const out = { a: m.a, b: m.b, section: m.section || defaultSection };
    if (m.material || defaultMaterial) out.material = m.material || defaultMaterial;
    if (m.group !== undefined) out.group = m.group;
    used[m.a] = 1; used[m.b] = 1;
    members.push(out);
  });
  let sup;
  if (supports === 'lowest') {
    const tol = supportTolerance === undefined ? Math.max(1e-3, 1e-3 * (zmax - zmin)) : supportTolerance;
    sup = [];
    for (let i = 0; i < n; i++) if (used[i] && nodes[i][2] <= zmin + tol) sup.push({ node: i, fix: supportType });
  } else if (Array.isArray(supports)) sup = supports;
  else throw new RangeError('supports must be "lowest" or an array of supports');
  return {
    nodes, members, supports: sup, loadCases: [{ name: 'selfWeight', selfWeight: {} }],
    meta: { kind: 'imported', name: frame.name, unusedNodes: Array.from(used).reduce((s, v) => s + (v ? 0 : 1), 0) },
  };
}
