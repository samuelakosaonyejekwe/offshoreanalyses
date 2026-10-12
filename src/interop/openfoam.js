// Hand-off to the finite-volume CFD toolbox OpenFOAM and back.
//
//   writeOpenfoamWaveTank     numerical wave tank around a fixed structure (interFoam, volume of fluid, wave
//                             generation and absorption boundary conditions)
//   writeOpenfoamCurrentCase  steady current past a fully submerged fixed body (simpleFoam, k-omega SST)
//   readOpenfoamForces        postProcessing/forces/<time>/force.dat (+ moment.dat) or the older forces.dat
//
// Target: OpenFOAM v2406 (openfoam.com, ESI line; v2312 uses the same dictionaries). The wave boundary
// conditions (waveVelocity, waveAlpha, constant/waveProperties with `waveModel`) belong to that line. The
// Foundation releases (v11, v12) differ in: constant/physicalProperties and momentumTransport instead of
// transportProperties and turbulenceProperties, their own wave library and waveProperties layout, and
// `foamRun -solver incompressibleVoF`; blockMeshDict, snappyHexMeshDict, decomposeParDict, setFieldsDict,
// constant/g and the forces function object are common to both lines.
//
// Coordinates: the app's (x, y, z up, still-water level z = 0). The tank bottom is at z = -waterDepth, waves
// travel towards +x, the structure stays where the model puts it.

import { worldCoords } from '../geometry/exporters.js';

const num = (v) => { if (!Number.isFinite(v)) throw new RangeError('openfoam: a non-finite number reached a dictionary'); return String(Number(v.toPrecision(10))); };
const vec = (v) => `(${num(v[0])} ${num(v[1])} ${num(v[2])})`;
const needPos = (v, name) => { if (!(v > 0) || !Number.isFinite(v)) throw new RangeError(`${name} must be > 0`); return v; };

function header(cls, object, location) {
  return ['FoamFile', '{', '    version     2.0;', '    format      ascii;', `    class       ${cls};`, ...(location ? [`    location    "${location}";`] : []), `    object      ${object};`, '}', ''].join('\n');
}
const dict = (object, location, body) => `${header('dictionary', object, location)}\n${body.trim()}\n`;

/** Field file: `patches` maps a patch name to the lines of its entry. */
function field(cls, object, dimensions, internal, patches) {
  const L = [header(cls, object, '0'), `dimensions      ${dimensions};`, '', `internalField   uniform ${internal};`, '', 'boundaryField', '{'];
  for (const [name, entry] of Object.entries(patches)) { L.push(`    ${name}`, '    {'); for (const e of entry) L.push(`        ${e}`); L.push('    }'); }
  L.push('}', '');
  return L.join('\n');
}

function surfaceInfo(surface) {
  if (!surface || surface.kind !== 'surface' || !surface.indices || surface.indices.length < 3) throw new RangeError("surface must be a 'surface' model with triangles");
  const P = worldCoords(surface);
  const lo = [Infinity, Infinity, Infinity]; const hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < P.length; i += 3) for (let k = 0; k < 3; k++) { if (P[i + k] < lo[k]) lo[k] = P[i + k]; if (P[i + k] > hi[k]) hi[k] = P[i + k]; }
  return { P, lo, hi, centre: [(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2] };
}

function triNormal(P, a, b, c) {
  const ux = P[b] - P[a]; const uy = P[b + 1] - P[a + 1]; const uz = P[b + 2] - P[a + 2];
  const vx = P[c] - P[a]; const vy = P[c + 1] - P[a + 1]; const vz = P[c + 2] - P[a + 2];
  const n = [uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx];
  const l = Math.hypot(n[0], n[1], n[2]) || 1;
  return [n[0] / l, n[1] / l, n[2] / l];
}

/** STL of the structure as one solid called `structure`: `{ text }` (ASCII) or `{ bytes }` (binary). */
function stlFile(P, indices, format) {
  const nT = indices.length / 3;
  if (format === 'ascii') {
    const L = ['solid structure'];
    for (let f = 0; f < nT; f++) {
      const a = 3 * indices[3 * f]; const b = 3 * indices[3 * f + 1]; const c = 3 * indices[3 * f + 2];
      const n = triNormal(P, a, b, c);
      L.push(` facet normal ${num(n[0])} ${num(n[1])} ${num(n[2])}`, '  outer loop');
      for (const q of [a, b, c]) L.push(`   vertex ${num(P[q])} ${num(P[q + 1])} ${num(P[q + 2])}`);
      L.push('  endloop', ' endfacet');
    }
    L.push('endsolid structure', '');
    return { text: L.join('\n') };
  }
  const bytes = new Uint8Array(84 + 50 * nT);
  bytes.set(new TextEncoder().encode('structure (binary STL, metres)'));
  const dv = new DataView(bytes.buffer);
  dv.setUint32(80, nT, true);
  for (let f = 0, o = 84; f < nT; f++, o += 50) {
    const a = 3 * indices[3 * f]; const b = 3 * indices[3 * f + 1]; const c = 3 * indices[3 * f + 2];
    const n = triNormal(P, a, b, c);
    for (let k = 0; k < 3; k++) dv.setFloat32(o + 4 * k, n[k], true);
    [a, b, c].forEach((q, v) => { for (let k = 0; k < 3; k++) dv.setFloat32(o + 12 + 12 * v + 4 * k, P[q + k], true); });
  }
  return { bytes };
}

function blockMeshDict(lo, hi, n, patches) {
  const v = [[lo[0], lo[1], lo[2]], [hi[0], lo[1], lo[2]], [hi[0], hi[1], lo[2]], [lo[0], hi[1], lo[2]], [lo[0], lo[1], hi[2]], [hi[0], lo[1], hi[2]], [hi[0], hi[1], hi[2]], [lo[0], hi[1], hi[2]]];
  const faces = { xMin: '(0 4 7 3)', xMax: '(1 2 6 5)', yMin: '(0 1 5 4)', yMax: '(3 7 6 2)', zMin: '(0 3 2 1)', zMax: '(4 5 6 7)' };
  const L = ['scale   1;', '', 'vertices', '(', ...v.map((p) => `    ${vec(p)}`), ');', '', 'blocks', '(', `    hex (0 1 2 3 4 5 6 7) (${n[0]} ${n[1]} ${n[2]}) simpleGrading (1 1 1)`, ');', '', 'edges', '(', ');', '', 'boundary', '('];
  for (const p of patches) { L.push(`    ${p.name}`, '    {', `        type ${p.type};`, '        faces', '        ('); for (const f of p.faces) L.push(`            ${faces[f]}`); L.push('        );', '    }'); }
  L.push(');', '', 'mergePatchPairs', '(', ');');
  return dict('blockMeshDict', 'system', L.join('\n'));
}

function snappyDict({ stlName, boxes, structureLevel, locationInMesh }) {
  const geo = [`    ${stlName}`, '    {', '        type triSurfaceMesh;', '        name structure;', '    }'];
  const reg = [];
  for (const b of boxes) {
    geo.push(`    ${b.name}`, '    {', '        type searchableBox;', `        min ${vec(b.min)};`, `        max ${vec(b.max)};`, '    }');
    reg.push(`        ${b.name}`, '        {', '            mode inside;', `            levels ((1e15 ${b.level}));`, '        }');
  }
  return dict('snappyHexMeshDict', 'system', `
castellatedMesh true;
snap            true;
addLayers       false;

geometry
{
${geo.join('\n')}
}

castellatedMeshControls
{
    maxLocalCells       2000000;
    maxGlobalCells      8000000;
    minRefinementCells  0;
    maxLoadUnbalance    0.10;
    nCellsBetweenLevels 3;
    features            ();
    refinementSurfaces
    {
        structure
        {
            level (${structureLevel} ${structureLevel});
            patchInfo
            {
                type wall;
            }
        }
    }
    resolveFeatureAngle 30;
    refinementRegions
    {
${reg.join('\n')}
    }
    locationInMesh ${vec(locationInMesh)};
    allowFreeStandingZoneFaces true;
}

snapControls
{
    nSmoothPatch    3;
    tolerance       2.0;
    nSolveIter      50;
    nRelaxIter      5;
    nFeatureSnapIter 10;
    implicitFeatureSnap true;
    explicitFeatureSnap false;
    multiRegionFeatureSnap false;
}

addLayersControls
{
    relativeSizes   true;
    layers
    {
    }
    expansionRatio  1.2;
    finalLayerThickness 0.5;
    minThickness    0.1;
    nGrow           0;
    featureAngle    60;
    nRelaxIter      5;
    nSmoothSurfaceNormals 1;
    nSmoothNormals  3;
    nSmoothThickness 10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3;
    minMedialAxisAngle 90;
    nBufferCellsNoExtrude 0;
    nLayerIter      50;
}

meshQualityControls
{
    maxNonOrtho     65;
    maxBoundarySkewness 20;
    maxInternalSkewness 4;
    maxConcave      80;
    minVol          1e-13;
    minTetQuality   -1e30;
    minArea         -1;
    minTwist        0.02;
    minDeterminant  0.001;
    minFaceWeight   0.05;
    minVolRatio     0.01;
    minTriangleTwist -1;
    nSmoothScale    4;
    errorReduction  0.75;
    relaxed
    {
        maxNonOrtho 75;
    }
}

mergeTolerance 1e-6;
`);
}

const decompose = (n) => dict('decomposeParDict', 'system', `numberOfSubdomains ${n};\n\nmethod          scotch;`);
// single-phase solvers work with kinematic pressure and need `rho rhoInf`; interFoam carries its own rho field
const forcesFunction = (rho, cofr, kinematic) => `    forces
    {
        type            forces;
        libs            (forces);
        patches         (structure);
        rhoInf          ${num(rho)};${kinematic ? '\n        rho             rhoInf;' : ''}
        log             on;
        writeControl    timeStep;
        writeInterval   1;
        CofR            ${vec(cofr)};
    }`;

const turbulenceDict = (model) => dict('turbulenceProperties', 'constant', model === 'kOmegaSST' ? `simulationType  RAS;

RAS
{
    RASModel        kOmegaSST;
    turbulence      on;
    printCoeffs     on;
}` : 'simulationType  laminar;');

/** A point in the fluid for snappyHexMesh: off the cell faces, away from the structure's bounding box. */
function pickLocation(lo, hi, cell, s) {
  const cand = [];
  for (const fx of [0.0371, 0.9629]) for (const fy of [0.0413, 0.9587]) for (const fz of [0.0457, 0.9543]) cand.push([lo[0] + fx * (hi[0] - lo[0]) + 0.137 * cell, lo[1] + fy * (hi[1] - lo[1]) + 0.119 * cell, lo[2] + fz * (hi[2] - lo[2]) + 0.123 * cell].map((v, k) => Math.min(Math.max(v, lo[k] + 0.31 * cell), hi[k] - 0.31 * cell)));
  for (const c of cand) if ([0, 1, 2].some((k) => c[k] < s.lo[k] - 0.5 * cell || c[k] > s.hi[k] + 0.5 * cell)) return c;
  throw new RangeError('domain: the structure fills the whole domain; no fluid point could be chosen for the mesher');
}

const WAVE_MODELS = { airy: 'StokesI', linear: 'StokesI', stokes1: 'StokesI', stokesI: 'StokesI', stokes2: 'StokesII', stokesII: 'StokesII', stokes5: 'StokesV', stokesV: 'StokesV', cnoidal: 'cnoidal' };

/**
 * Numerical wave tank around a fixed structure for interFoam (OpenFOAM v2406 / v2312): regular waves are
 * generated at the inlet (x-min) and absorbed at the outlet (x-max) by the wave boundary conditions; the forces
 * and moments on the structure are logged every time step.
 *
 * Files: constant/triSurface/<name>.stl, system/{blockMeshDict, snappyHexMeshDict, controlDict, fvSchemes,
 * fvSolution, decomposeParDict, setFieldsDict}, constant/{transportProperties, turbulenceProperties, g,
 * waveProperties} (+ the other turbulence variant as turbulenceProperties.<variant>), 0/{alpha.water,
 * alpha.water.orig, U, p_rgh, k, omega, nut}, Allrun, Allclean.
 * Patches: inlet, outlet, sides (slip), bottom (wall), atmosphere, structure (wall, from the STL).
 * Reference: OpenFOAM v2406 User Guide (blockMesh, snappyHexMesh, interFoam); Higuera, Lara & Losada, "Realistic
 * wave generation and active wave absorption for Navier–Stokes models", Coastal Eng. 71 (2013) (the wave
 * boundary conditions); Hirt & Nichols, J. Comput. Phys. 39 (1981) (volume of fluid).
 * Validity: fixed (non-moving) structure; long-crested regular waves along +x; the outlet absorbs with
 * shallow-water theory, so reflections grow for short waves in deep tanks (use a longer tank and stop the run
 * before the reflection returns). No current: these wave boundary conditions do not superpose one (a non-zero
 * `current` is reported in the notes and ignored). Resolve the wave height with about 10 cells and the wave
 * length with about 80 before trusting forces.
 * @param {object} o
 * @param {object} o.surface structure surface (CONTRACT §4 `surface`), m, z up, still water at z = 0
 * @param {{length:number, width:number, height:number, waterDepth:number, xMin?:number, yMin?:number}} o.domain
 *   tank length (x), width (y), total height (z) and still-water depth, m; xMin (default: one third of the
 *   length up-wave of the structure centre) and yMin (default: centred on the structure), m
 * @param {{H:number, T:number, theory?:string, rampTime?:number}} o.wave height m, period s, theory 'airy' |
 *   'stokes2' | 'stokes5' | 'cnoidal' (default 'stokes2'), ramp time s (default one period)
 * @param {number} [o.current=0] m/s — not applied, see above
 * @param {number} o.cellSize background cell size, m
 * @param {{structure?:number, box?:number, freeSurface?:number, boxMargin?:number}} [o.refinement] refinement
 *   levels (each halves the cell): on the structure (default 2), in a box around it (1) and in the band
 *   |z| ≤ H around the free surface (1); boxMargin: clearance of the box around the structure, m (default 4 cells)
 * @param {number} o.endTime simulated time, s
 * @param {'laminar'|'kOmegaSST'} [o.turbulence='laminar'] active turbulence model
 * @param {number} [o.processors=4] subdomains for the parallel run
 * @param {'binary'|'ascii'} [o.stl='binary'] STL encoding
 * @param {number[]} [o.cofr] centre about which moments are reported, m (default: the origin)
 * @param {number} [o.writeInterval] field output interval, s (default T / 10)
 * @param {{rhoWater?:number, nuWater?:number, rhoAir?:number, nuAir?:number, sigma?:number, g?:number}} [o.fluid]
 *   kg/m³, m²/s, kg/m³, m²/s, N/m, m/s² (defaults 1025, 1.05e-6, 1.2, 1.48e-5, 0.07, 9.80665)
 * @param {string} [o.name] STL file name (default: the surface name)
 * @returns {{files:{path:string, text?:string, bytes?:Uint8Array}[], notes:string[], mesh:{background:number[],
 *   backgroundCells:number, bounds:{min:number[], max:number[]}}}} case folder; background cell counts
 */
export function writeOpenfoamWaveTank({ surface, domain, wave, current = 0, cellSize, refinement = {}, endTime, turbulence = 'laminar', processors = 4, stl = 'binary', cofr = [0, 0, 0], writeInterval, fluid = {}, name } = {}) {
  const s = surfaceInfo(surface);
  if (!domain) throw new RangeError('domain is required');
  needPos(domain.length, 'domain.length'); needPos(domain.width, 'domain.width'); needPos(domain.height, 'domain.height'); needPos(domain.waterDepth, 'domain.waterDepth');
  if (!(domain.height > domain.waterDepth)) throw new RangeError('domain.height must exceed domain.waterDepth (air is needed above the water)');
  if (!wave) throw new RangeError('wave is required');
  needPos(wave.H, 'wave.H'); needPos(wave.T, 'wave.T'); needPos(cellSize, 'cellSize'); needPos(endTime, 'endTime');
  const model = WAVE_MODELS[wave.theory === undefined ? 'stokes2' : wave.theory];
  if (!model) throw new RangeError(`wave.theory must be one of ${Object.keys(WAVE_MODELS).join(', ')}`);
  if (!(wave.H < 0.78 * domain.waterDepth)) throw new RangeError('wave.H exceeds the depth-limited breaking height 0.78 × waterDepth');
  if (!(domain.height - domain.waterDepth > wave.H)) throw new RangeError('domain.height: leave more than one wave height of air above still water');
  if (!['laminar', 'kOmegaSST'].includes(turbulence)) throw new RangeError("turbulence must be 'laminar' or 'kOmegaSST'");
  if (!['binary', 'ascii'].includes(stl)) throw new RangeError("stl must be 'binary' or 'ascii'");
  if (!Number.isInteger(processors) || processors < 1) throw new RangeError('processors must be a positive integer');
  const lvl = { structure: refinement.structure ?? 2, box: refinement.box ?? 1, freeSurface: refinement.freeSurface ?? 1 };
  for (const k of Object.keys(lvl)) if (!Number.isInteger(lvl[k]) || lvl[k] < 0 || lvl[k] > 6) throw new RangeError(`refinement.${k} must be an integer level 0…6`);
  const fl = { rhoWater: 1025, nuWater: 1.05e-6, rhoAir: 1.2, nuAir: 1.48e-5, sigma: 0.07, g: 9.80665, ...fluid };
  const nm = String(name || surface.name || 'structure').replace(/[^A-Za-z0-9_-]+/g, '_') || 'structure';
  const xMin = domain.xMin ?? s.centre[0] - domain.length / 3;
  const yMin = domain.yMin ?? s.centre[1] - domain.width / 2;
  const lo = [xMin, yMin, -domain.waterDepth]; const hi = [xMin + domain.length, yMin + domain.width, domain.height - domain.waterDepth];
  const n = [0, 1, 2].map((k) => Math.max(1, Math.round((hi[k] - lo[k]) / cellSize)));
  const notes = [];
  if (s.lo[0] <= lo[0] || s.hi[0] >= hi[0] || s.lo[1] <= lo[1] || s.hi[1] >= hi[1]) throw new RangeError('domain: the structure reaches the inlet, outlet or side walls; enlarge the tank or move it (xMin, yMin)');
  if (s.lo[2] < lo[2]) notes.push('The structure extends below the tank bottom; the part below z = -waterDepth is cut off by the mesher.');
  if (s.hi[2] > hi[2]) notes.push('The structure extends above the top of the tank; the part above is cut off by the mesher.');
  if (current) notes.push(`current = ${num(current)} m/s is NOT applied: the wave generation boundary conditions of OpenFOAM v2406 do not superpose a current. Use writeOpenfoamCurrentCase for current loads.`);
  const margin = refinement.boxMargin ?? 4 * cellSize;
  const boxes = [];
  if (lvl.box > 0) boxes.push({ name: 'refinementBox', level: lvl.box, min: s.lo.map((v, k) => Math.max(v - margin, lo[k])), max: s.hi.map((v, k) => Math.min(v + margin, hi[k])) });
  if (lvl.freeSurface > 0) boxes.push({ name: 'freeSurface', level: lvl.freeSurface, min: [lo[0], lo[1], Math.max(-wave.H, lo[2])], max: [hi[0], hi[1], Math.min(wave.H, hi[2])] });
  const fine = cellSize / 2 ** Math.max(lvl.freeSurface, 0);
  const cellsPerH = wave.H / fine;
  if (cellsPerH < 8) notes.push(`Only ${num(Number(cellsPerH.toPrecision(3)))} cells per wave height at the free surface (about 10 are needed for forces): reduce cellSize or raise refinement.freeSurface.`);
  const wInt = writeInterval ?? wave.T / 10;
  const ramp = wave.rampTime ?? wave.T;
  const wallBC = (kind) => ({ U: ['type            noSlip;'], k: ['type            kqRWallFunction;', 'value           uniform 1e-06;'], omega: ['type            omegaWallFunction;', 'value           uniform 1;'], nut: ['type            nutkWallFunction;', 'value           uniform 0;'] })[kind];
  const alphaPatches = {
    inlet: ['type            waveAlpha;', 'value           uniform 0;'], outlet: ['type            zeroGradient;'], sides: ['type            zeroGradient;'], bottom: ['type            zeroGradient;'],
    atmosphere: ['type            inletOutlet;', 'inletValue      uniform 0;', 'value           uniform 0;'], structure: ['type            zeroGradient;'],
  };
  const flux = ['type            fixedFluxPressure;', 'value           uniform 0;'];
  const scalarTurb = (v, wall) => ({ inlet: ['type            fixedValue;', `value           uniform ${v};`], outlet: ['type            inletOutlet;', `inletValue      uniform ${v};`, `value           uniform ${v};`], sides: ['type            zeroGradient;'], bottom: wall, atmosphere: ['type            inletOutlet;', `inletValue      uniform ${v};`, `value           uniform ${v};`], structure: wall });
  const files = [];
  const add = (path, text) => files.push({ path, text });
  files.push({ path: `constant/triSurface/${nm}.stl`, ...stlFile(s.P, surface.indices, stl) });
  add('system/blockMeshDict', blockMeshDict(lo, hi, n, [{ name: 'inlet', type: 'patch', faces: ['xMin'] }, { name: 'outlet', type: 'patch', faces: ['xMax'] }, { name: 'sides', type: 'patch', faces: ['yMin', 'yMax'] }, { name: 'bottom', type: 'wall', faces: ['zMin'] }, { name: 'atmosphere', type: 'patch', faces: ['zMax'] }]));
  add('system/snappyHexMeshDict', snappyDict({ stlName: `${nm}.stl`, boxes, structureLevel: lvl.structure, locationInMesh: pickLocation(lo, hi, cellSize, s) }));
  add('system/controlDict', dict('controlDict', 'system', `
application     interFoam;

startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         ${num(endTime)};
deltaT          ${num(Math.min(wave.T / 400, endTime / 10))};

writeControl    adjustable;
writeInterval   ${num(wInt)};
purgeWrite      0;
writeFormat     ascii;
writePrecision  6;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable yes;

adjustTimeStep  yes;
maxCo           0.5;
maxAlphaCo      0.5;
maxDeltaT       ${num(wave.T / 50)};

functions
{
${forcesFunction(fl.rhoWater, cofr, false)}
}`));
  add('system/fvSchemes', dict('fvSchemes', 'system', `
ddtSchemes
{
    default         Euler;
}

gradSchemes
{
    default         Gauss linear;
    limitedGrad     cellLimited Gauss linear 1;
}

divSchemes
{
    div(rhoPhi,U)   Gauss linearUpwind grad(U);
    div(phi,alpha)  Gauss vanLeer;
    div(phirb,alpha) Gauss linear;
    div(phi,k)      Gauss linearUpwind limitedGrad;
    div(phi,omega)  Gauss linearUpwind limitedGrad;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}

wallDist
{
    method          meshWave;
}`));
  add('system/fvSolution', dict('fvSolution', 'system', `
solvers
{
    "alpha.water.*"
    {
        nAlphaCorr      1;
        nAlphaSubCycles 2;
        cAlpha          1;
    }

    "pcorr.*"
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-6;
        relTol          0;
    }

    p_rgh
    {
        solver          GAMG;
        smoother        DIC;
        tolerance       1e-7;
        relTol          0.05;
    }

    p_rghFinal
    {
        solver          GAMG;
        smoother        DIC;
        tolerance       1e-8;
        relTol          0;
    }

    "(U|k|omega)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-7;
        relTol          0.1;
    }

    "(U|k|omega)Final"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-7;
        relTol          0;
    }
}

PIMPLE
{
    momentumPredictor no;
    nOuterCorrectors 1;
    nCorrectors     3;
    nNonOrthogonalCorrectors 1;
}

relaxationFactors
{
    equations
    {
        ".*"            1;
    }
}`));
  add('system/decomposeParDict', decompose(processors));
  add('system/setFieldsDict', dict('setFieldsDict', 'system', `
defaultFieldValues
(
    volScalarFieldValue alpha.water 0
);

regions
(
    boxToCell
    {
        box ${vec([lo[0] - 1, lo[1] - 1, lo[2] - 1])} ${vec([hi[0] + 1, hi[1] + 1, 0])};
        fieldValues
        (
            volScalarFieldValue alpha.water 1
        );
    }
);`));
  add('constant/transportProperties', dict('transportProperties', 'constant', `
phases          (water air);

water
{
    transportModel  Newtonian;
    nu              ${num(fl.nuWater)};
    rho             ${num(fl.rhoWater)};
}

air
{
    transportModel  Newtonian;
    nu              ${num(fl.nuAir)};
    rho             ${num(fl.rhoAir)};
}

sigma           ${num(fl.sigma)};`));
  const other = turbulence === 'laminar' ? 'kOmegaSST' : 'laminar';
  add('constant/turbulenceProperties', turbulenceDict(turbulence));
  add(`constant/turbulenceProperties.${other}`, turbulenceDict(other));
  add('constant/g', `${header('uniformDimensionedVectorField', 'g', 'constant')}\ndimensions      [0 1 -2 0 0 0 0];\nvalue           (0 0 ${num(-fl.g)});\n`);
  add('constant/waveProperties', dict('waveProperties', 'constant', `
inlet
{
    alpha           alpha.water;
    waveModel       ${model};
    nPaddle         1;
    waveHeight      ${num(wave.H)};
    waveAngle       0.0;
    rampTime        ${num(ramp)};
    activeAbsorption yes;
    wavePeriod      ${num(wave.T)};
}

outlet
{
    alpha           alpha.water;
    waveModel       shallowWaterAbsorption;
    nPaddle         1;
}`));
  const alphaText = field('volScalarField', 'alpha.water', '[0 0 0 0 0 0 0]', '0', alphaPatches);
  add('0/alpha.water', alphaText);
  add('0/alpha.water.orig', alphaText);
  add('0/U', field('volVectorField', 'U', '[0 1 -1 0 0 0 0]', '(0 0 0)', {
    inlet: ['type            waveVelocity;', 'value           uniform (0 0 0);'], outlet: ['type            waveVelocity;', 'value           uniform (0 0 0);'], sides: ['type            slip;'],
    bottom: wallBC('U'), atmosphere: ['type            pressureInletOutletVelocity;', 'value           uniform (0 0 0);'], structure: wallBC('U'),
  }));
  add('0/p_rgh', field('volScalarField', 'p_rgh', '[1 -1 -2 0 0 0 0]', '0', { inlet: flux, outlet: flux, sides: flux, bottom: flux, atmosphere: ['type            totalPressure;', 'p0              uniform 0;'], structure: flux }));
  add('0/k', field('volScalarField', 'k', '[0 2 -2 0 0 0 0]', '1e-06', scalarTurb('1e-06', wallBC('k'))));
  add('0/omega', field('volScalarField', 'omega', '[0 0 -1 0 0 0 0]', '1', scalarTurb('1', wallBC('omega'))));
  add('0/nut', field('volScalarField', 'nut', '[0 2 -1 0 0 0 0]', '0', { inlet: ['type            calculated;', 'value           uniform 0;'], outlet: ['type            calculated;', 'value           uniform 0;'], sides: ['type            calculated;', 'value           uniform 0;'], bottom: wallBC('nut'), atmosphere: ['type            calculated;', 'value           uniform 0;'], structure: wallBC('nut') }));
  add('Allrun', `#!/bin/sh
cd "\${0%/*}" || exit 1
. "\${WM_PROJECT_DIR:?}/bin/tools/RunFunctions"

runApplication blockMesh
runApplication snappyHexMesh -overwrite
runApplication checkMesh
cp 0/alpha.water.orig 0/alpha.water
runApplication setFields
${processors > 1 ? `runApplication decomposePar
runParallel $(getApplication)
runApplication reconstructPar` : 'runApplication $(getApplication)'}
`);
  add('Allclean', `#!/bin/sh
cd "\${0%/*}" || exit 1
. "\${WM_PROJECT_DIR:?}/bin/tools/CleanFunctions"

cleanTimeDirectories
rm -rf constant/polyMesh constant/extendedFeatureEdgeMesh processor* postProcessing log.* dynamicCode
cp 0/alpha.water.orig 0/alpha.water
`);
  notes.push(`Background mesh ${n[0]} × ${n[1]} × ${n[2]} = ${n[0] * n[1] * n[2]} cells of ${num(cellSize)} m before refinement; wave model ${model}, H = ${num(wave.H)} m, T = ${num(wave.T)} s, still-water depth ${num(domain.waterDepth)} m.`);
  notes.push('Run:  chmod +x Allrun Allclean ; ./Allrun   in an OpenFOAM v2406 (or v2312) shell. Forces: postProcessing/forces/0/force.dat and moment.dat, read with readOpenfoamForces.');
  notes.push(`Turbulence: ${turbulence}; to switch, copy constant/turbulenceProperties.${other} over constant/turbulenceProperties (the k, omega and nut fields are already there).`);
  notes.push('Mesh quality (checkMesh) and time-step/mesh convergence remain the user\'s responsibility; the case is a starting point, not a validated set-up.');
  return { files, notes, mesh: { background: n, backgroundCells: n[0] * n[1] * n[2], bounds: { min: lo, max: hi } } };
}

/**
 * Steady uniform current past a fixed, fully submerged body for simpleFoam (OpenFOAM v2406 / v2312, k-omega SST
 * or laminar): inlet at x-min with the current speed, pressure outlet at x-max, slip on the other four faces
 * (no free surface is modelled), no-slip wall with wall functions on the structure, forces logged every iteration.
 * Reference: OpenFOAM v2406 User Guide (simpleFoam); Menter, Kuntz & Langtry, Turbulence, Heat and Mass
 * Transfer 4 (2003) (k-omega SST).
 * Validity: single-phase, incompressible, steady RANS; bluff bodies shed vortices and need a transient solver
 * for anything but a mean-load estimate; blockage (body frontal area / domain cross-section) should stay below
 * about 5 %.
 * @param {object} o
 * @param {object} o.surface structure surface, m
 * @param {{length:number, width:number, height:number, xMin?:number, yMin?:number, zMin?:number}} o.domain box
 *   dimensions, m; default position: body one third of the length from the inlet, centred in y and z
 * @param {number} o.current current speed along +x, m/s (> 0)
 * @param {number} o.cellSize background cell size, m
 * @param {{structure?:number, box?:number, boxMargin?:number}} [o.refinement] levels on the structure (default 2)
 *   and in a box around it (1); box clearance, m (default 4 cells)
 * @param {number} [o.iterations=1000] SIMPLE iterations
 * @param {'kOmegaSST'|'laminar'} [o.turbulence='kOmegaSST']
 * @param {number} [o.turbulenceIntensity=0.05] inlet turbulence intensity (-)
 * @param {{rho?:number, nu?:number}} [o.fluid] kg/m³, m²/s (defaults 1025, 1.05e-6)
 * @param {number} [o.processors=4]
 * @param {'binary'|'ascii'} [o.stl='binary']
 * @param {number[]} [o.cofr] moment reference point, m (default: centre of the body's bounding box)
 * @param {string} [o.name] STL file name
 * @returns {{files:{path:string, text?:string, bytes?:Uint8Array}[], notes:string[], mesh:{background:number[],
 *   backgroundCells:number, bounds:{min:number[], max:number[]}}}}
 */
export function writeOpenfoamCurrentCase({ surface, domain, current, cellSize, refinement = {}, iterations = 1000, turbulence = 'kOmegaSST', turbulenceIntensity = 0.05, fluid = {}, processors = 4, stl = 'binary', cofr, name } = {}) {
  const s = surfaceInfo(surface);
  if (!domain) throw new RangeError('domain is required');
  needPos(domain.length, 'domain.length'); needPos(domain.width, 'domain.width'); needPos(domain.height, 'domain.height');
  needPos(current, 'current'); needPos(cellSize, 'cellSize'); needPos(turbulenceIntensity, 'turbulenceIntensity');
  if (!Number.isInteger(iterations) || iterations < 1) throw new RangeError('iterations must be a positive integer');
  if (!['laminar', 'kOmegaSST'].includes(turbulence)) throw new RangeError("turbulence must be 'laminar' or 'kOmegaSST'");
  if (!['binary', 'ascii'].includes(stl)) throw new RangeError("stl must be 'binary' or 'ascii'");
  if (!Number.isInteger(processors) || processors < 1) throw new RangeError('processors must be a positive integer');
  const lvl = { structure: refinement.structure ?? 2, box: refinement.box ?? 1 };
  for (const k of Object.keys(lvl)) if (!Number.isInteger(lvl[k]) || lvl[k] < 0 || lvl[k] > 6) throw new RangeError(`refinement.${k} must be an integer level 0…6`);
  const fl = { rho: 1025, nu: 1.05e-6, ...fluid };
  const nm = String(name || surface.name || 'structure').replace(/[^A-Za-z0-9_-]+/g, '_') || 'structure';
  const lo = [domain.xMin ?? s.centre[0] - domain.length / 3, domain.yMin ?? s.centre[1] - domain.width / 2, domain.zMin ?? s.centre[2] - domain.height / 2];
  const hi = [lo[0] + domain.length, lo[1] + domain.width, lo[2] + domain.height];
  for (let k = 0; k < 3; k++) if (s.lo[k] <= lo[k] || s.hi[k] >= hi[k]) throw new RangeError('domain: the body must lie strictly inside the box (enlarge it or set xMin, yMin, zMin)');
  const n = [0, 1, 2].map((k) => Math.max(1, Math.round((hi[k] - lo[k]) / cellSize)));
  const margin = refinement.boxMargin ?? 4 * cellSize;
  const boxes = lvl.box > 0 ? [{ name: 'refinementBox', level: lvl.box, min: s.lo.map((v, k) => Math.max(v - margin, lo[k])), max: s.hi.map((v, k) => Math.min(v + margin, hi[k])) }] : [];
  const size = Math.max(s.hi[1] - s.lo[1], s.hi[2] - s.lo[2]);
  const kIn = 1.5 * (turbulenceIntensity * current) ** 2;
  // omega = sqrt(k) / (C_mu^(1/4) l), mixing length l = 7 % of the body's cross-flow size
  const omIn = Math.sqrt(kIn) / (0.09 ** 0.25 * 0.07 * size);
  const centre = cofr || s.centre;
  const notes = [];
  const blockage = ((s.hi[1] - s.lo[1]) * (s.hi[2] - s.lo[2])) / (domain.width * domain.height);
  if (blockage > 0.05) notes.push(`Blockage ${num(Number((100 * blockage).toPrecision(3)))} % (bounding-box frontal area over the domain cross-section) exceeds 5 %: the forces will be overestimated; widen the domain.`);
  const slip = ['type            slip;'];
  const files = [];
  const add = (path, text) => files.push({ path, text });
  files.push({ path: `constant/triSurface/${nm}.stl`, ...stlFile(s.P, surface.indices, stl) });
  add('system/blockMeshDict', blockMeshDict(lo, hi, n, [{ name: 'inlet', type: 'patch', faces: ['xMin'] }, { name: 'outlet', type: 'patch', faces: ['xMax'] }, { name: 'sides', type: 'patch', faces: ['yMin', 'yMax', 'zMin', 'zMax'] }]));
  add('system/snappyHexMeshDict', snappyDict({ stlName: `${nm}.stl`, boxes, structureLevel: lvl.structure, locationInMesh: pickLocation(lo, hi, cellSize, s) }));
  add('system/controlDict', dict('controlDict', 'system', `
application     simpleFoam;

startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         ${iterations};
deltaT          1;

writeControl    timeStep;
writeInterval   ${Math.max(1, Math.round(iterations / 4))};
purgeWrite      2;
writeFormat     ascii;
writePrecision  6;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable yes;

functions
{
${forcesFunction(fl.rho, centre, true)}
}`));
  add('system/fvSchemes', dict('fvSchemes', 'system', `
ddtSchemes
{
    default         steadyState;
}

gradSchemes
{
    default         Gauss linear;
    limitedGrad     cellLimited Gauss linear 1;
}

divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind limitedGrad;
    div(phi,k)      bounded Gauss upwind;
    div(phi,omega)  bounded Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}

wallDist
{
    method          meshWave;
}`));
  add('system/fvSolution', dict('fvSolution', 'system', `
solvers
{
    p
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-7;
        relTol          0.05;
    }

    "(U|k|omega)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0.1;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 1;
    consistent      yes;
    residualControl
    {
        p               1e-4;
        U               1e-5;
        "(k|omega)"     1e-5;
    }
}

relaxationFactors
{
    equations
    {
        U               0.9;
        ".*"            0.8;
    }
}`));
  add('system/decomposeParDict', decompose(processors));
  add('constant/transportProperties', dict('transportProperties', 'constant', `transportModel  Newtonian;\n\nnu              ${num(fl.nu)};`));
  const other = turbulence === 'laminar' ? 'kOmegaSST' : 'laminar';
  add('constant/turbulenceProperties', turbulenceDict(turbulence));
  add(`constant/turbulenceProperties.${other}`, turbulenceDict(other));
  const U = `(${num(current)} 0 0)`;
  const io = (v) => ['type            inletOutlet;', `inletValue      uniform ${v};`, `value           uniform ${v};`];
  add('0/U', field('volVectorField', 'U', '[0 1 -1 0 0 0 0]', U, { inlet: ['type            fixedValue;', `value           uniform ${U};`], outlet: io(U), sides: slip, structure: ['type            noSlip;'] }));
  add('0/p', field('volScalarField', 'p', '[0 2 -2 0 0 0 0]', '0', { inlet: ['type            zeroGradient;'], outlet: ['type            fixedValue;', 'value           uniform 0;'], sides: ['type            zeroGradient;'], structure: ['type            zeroGradient;'] }));
  add('0/k', field('volScalarField', 'k', '[0 2 -2 0 0 0 0]', num(kIn), { inlet: ['type            fixedValue;', `value           uniform ${num(kIn)};`], outlet: io(num(kIn)), sides: ['type            zeroGradient;'], structure: ['type            kqRWallFunction;', `value           uniform ${num(kIn)};`] }));
  add('0/omega', field('volScalarField', 'omega', '[0 0 -1 0 0 0 0]', num(omIn), { inlet: ['type            fixedValue;', `value           uniform ${num(omIn)};`], outlet: io(num(omIn)), sides: ['type            zeroGradient;'], structure: ['type            omegaWallFunction;', `value           uniform ${num(omIn)};`] }));
  add('0/nut', field('volScalarField', 'nut', '[0 2 -1 0 0 0 0]', '0', { inlet: ['type            calculated;', 'value           uniform 0;'], outlet: ['type            calculated;', 'value           uniform 0;'], sides: ['type            calculated;', 'value           uniform 0;'], structure: ['type            nutkWallFunction;', 'value           uniform 0;'] }));
  add('Allrun', `#!/bin/sh
cd "\${0%/*}" || exit 1
. "\${WM_PROJECT_DIR:?}/bin/tools/RunFunctions"

runApplication blockMesh
runApplication snappyHexMesh -overwrite
runApplication checkMesh
${processors > 1 ? `runApplication decomposePar
runParallel $(getApplication)
runApplication reconstructPar -latestTime` : 'runApplication $(getApplication)'}
`);
  add('Allclean', `#!/bin/sh
cd "\${0%/*}" || exit 1
. "\${WM_PROJECT_DIR:?}/bin/tools/CleanFunctions"

cleanTimeDirectories
rm -rf constant/polyMesh constant/extendedFeatureEdgeMesh processor* postProcessing log.* dynamicCode
`);
  notes.push(`Background mesh ${n[0]} × ${n[1]} × ${n[2]} = ${n[0] * n[1] * n[2]} cells of ${num(cellSize)} m before refinement; U = ${num(current)} m/s, k = ${num(Number(kIn.toPrecision(4)))} m²/s², omega = ${num(Number(omIn.toPrecision(4)))} 1/s at the inlet.`);
  notes.push('Run:  chmod +x Allrun Allclean ; ./Allrun   in an OpenFOAM v2406 (or v2312) shell. Forces: postProcessing/forces/0/force.dat and moment.dat (N, N·m; average the last iterations).');
  notes.push('No free surface: valid for bodies submerged by several diameters. Wall functions need y+ between about 30 and 300 on the structure; check it with the yPlus function object.');
  return { files, notes, mesh: { background: n, backgroundCells: n[0] * n[1] * n[2], bounds: { min: lo, max: hi } } };
}

/**
 * Parses the output of the OpenFOAM `forces` function object into time series.
 * Two layouts are recognised:
 *  - bracketed (OpenFOAM ≤ v2012 and the Foundation releases, `forces.dat`):
 *      `time ((px py pz) (vx vy vz) [(porous)]) ((mpx mpy mpz) (mvx mvy mvz) [(porous)])`
 *  - columns (OpenFOAM v2106 and later, `force.dat` and `moment.dat` as separate files):
 *      `time total_x total_y total_z pressure_x … viscous_x … [porous_x …]`
 * Restarted runs that repeat a time keep the later record.
 * Reference: OpenFOAM v2406 documentation, function object `forces` (output files); earlier format from the
 * v2012 / v7 sources (forces.C, writeFileHeader).
 * @param {string} text content of force.dat or forces.dat
 * @param {{momentText?:string}} [opts] content of moment.dat (column layout) to add the moments
 * @returns {{format:'bracketed'|'columns', cofr:number[]|null, t:number[], Fx:number[], Fy:number[], Fz:number[],
 *   Mx:number[]|null, My:number[]|null, Mz:number[]|null, pressure:{Fx:number[], Fy:number[], Fz:number[]},
 *   viscous:{Fx:number[], Fy:number[], Fz:number[]}, n:number}} time s; total force N and moment N·m about
 *   `cofr` (null when moments are not available); pressure and viscous parts of the force, N
 */
export function readOpenfoamForces(text, { momentText } = {}) {
  const parse = (src) => {
    const rows = []; let cofr = null; let bracket = false;
    for (const line of String(src).split(/\r?\n/)) {
      const t = line.trim();
      if (!t) continue;
      if (t[0] === '#') { const m = /CofR\s*:\s*\(([^)]*)\)/.exec(t); if (m) cofr = m[1].trim().split(/\s+/).map(Number); continue; }
      if (t.includes('(')) bracket = true;
      const v = t.replace(/[()]/g, ' ').trim().split(/\s+/).map(Number);
      if (v.length >= 4 && v.every((x) => Number.isFinite(x))) rows.push(v);
    }
    // later records win when a restart repeats times
    const byTime = new Map();
    for (const r of rows) byTime.set(r[0], r);
    return { rows: Array.from(byTime.values()).sort((a, b) => a[0] - b[0]), cofr, bracket };
  };
  const f = parse(text);
  if (!f.rows.length) throw new RangeError('forces: no data rows found (expected the force.dat or forces.dat written by the forces function object)');
  const col = (rows, k) => rows.map((r) => r[k]);
  const out = { format: f.bracket ? 'bracketed' : 'columns', cofr: f.cofr, t: col(f.rows, 0), n: f.rows.length, Mx: null, My: null, Mz: null };
  if (f.bracket) {
    const w = f.rows[0].length;
    // 13 columns: (pressure, viscous) × (force, moment); 19: with porous parts
    const groups = w === 13 ? 2 : w === 19 ? 3 : 0;
    if (!groups) throw new RangeError(`forces: bracketed record with ${w} numbers is not a known layout (expected 13 or 19)`);
    const g = 3 * groups;
    const tot = (k0) => f.rows.map((r) => { let s = 0; for (let q = 0; q < groups; q++) s += r[k0 + 3 * q]; return s; });
    out.Fx = tot(1); out.Fy = tot(2); out.Fz = tot(3);
    out.Mx = tot(1 + g); out.My = tot(2 + g); out.Mz = tot(3 + g);
    out.pressure = { Fx: col(f.rows, 1), Fy: col(f.rows, 2), Fz: col(f.rows, 3) };
    out.viscous = { Fx: col(f.rows, 4), Fy: col(f.rows, 5), Fz: col(f.rows, 6) };
    return out;
  }
  if (f.rows[0].length < 10) throw new RangeError('forces: column layout needs time, total, pressure and viscous components (10 columns)');
  out.Fx = col(f.rows, 1); out.Fy = col(f.rows, 2); out.Fz = col(f.rows, 3);
  out.pressure = { Fx: col(f.rows, 4), Fy: col(f.rows, 5), Fz: col(f.rows, 6) };
  out.viscous = { Fx: col(f.rows, 7), Fy: col(f.rows, 8), Fz: col(f.rows, 9) };
  if (momentText) {
    const m = parse(momentText);
    if (m.bracket) throw new RangeError('forces: momentText must be a column-layout moment.dat');
    const byTime = new Map(m.rows.map((r) => [r[0], r]));
    const pick = (k) => out.t.map((t) => { const r = byTime.get(t); return r ? r[k] : NaN; });
    out.Mx = pick(1); out.My = pick(2); out.Mz = pick(3);
    if (!out.cofr) out.cofr = m.cofr;
  }
  return out;
}
