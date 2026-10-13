// Geometry format registry: the compatibility matrix shown to users, format detection and import dispatch.
//
// Each entry of FORMATS says what really happens to a file of that format:
//   status 'native'  — read here, in plain JavaScript (`via` names the importer module)
//   status 'adapter' — needs a WebAssembly reader in src/geometry/adapters (`via: 'adapter:<name>'`); when the
//                      adapter is not part of the build, importFile() says so and names what to export instead.
//                      DWG is the one adapter whose library is not shipped at all: it is downloaded at the
//                      user's request (opts.consent) and runs in a worker of its own.
//                      Where `nativeSubset` is set, that part of the format is read without the adapter.
//   status 'convert' — not read at all; importFile() throws with the exact export to make in the source program.
// `suites` lists the suites (1–16) whose specification names the format; an empty list means the format is an
// extra that no suite asks for.

import { detectFormat, detectFileSet, EXTENSION_MAP } from './detect.js';
import { zipEntries, readZip, gunzip } from './zip.js';
import { extOf, finish } from './importers/common.js';

const ALL = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];
const pick = (promise, fn, extra) => promise.then((m) => ({ parse: extra ? (i, o) => m[fn](i, { ...o, ...extra }) : m[fn] }));

const ADAPTER_LOADERS = {
  occt: () => import('./adapters/occt.js'),
  'web-ifc': () => import('./adapters/web-ifc.js'),
  hdf5: () => import('./adapters/hdf5.js'),
  laz: () => import('./adapters/laz.js'),
  sqlite: () => import('./adapters/sqlite.js'),
};
// Readers that are not part of the build: their library is fetched at run time, at the user's request, by a worker
// of its own (see adapters/dwg.js). They are resolved like the others but are not listed by adapterStatus().
const REMOTE_LOADERS = {
  dwg: () => import('./adapters/dwg.js'),
};
const loaderOf = (name) => ADAPTER_LOADERS[name] || REMOTE_LOADERS[name];
const adapter = (name) => () => loaderOf(name)();

/** @type {ReadonlyArray<object>} */
export const FORMATS = Object.freeze([
  // ---- triangulated surfaces and scans ----
  { id: 'stl', name: 'STL (ASCII and binary)', extensions: ['stl'], kinds: ['surface'], category: 'mesh', status: 'native', via: 'importers/surfaces.js', suites: ALL,
    notes: 'Facet corners are welded into an indexed mesh; solids of an ASCII file become groups. No unit in the format: metres assumed unless stated.', load: () => pick(import('./importers/surfaces.js'), 'parseStl') },
  { id: 'obj', name: 'Wavefront OBJ', extensions: ['obj'], kinds: ['surface', 'polyline', 'points'], category: 'mesh', status: 'native', via: 'importers/surfaces.js', suites: ALL,
    notes: 'Faces (n-gons triangulated), lines and points; o / g / usemtl groups. Normals, textures and free-form curves are not kept.', load: () => pick(import('./importers/surfaces.js'), 'parseObj') },
  { id: 'ply', name: 'Stanford PLY (ASCII and binary)', extensions: ['ply'], kinds: ['surface', 'points'], category: 'mesh', status: 'native', via: 'importers/surfaces.js', suites: ALL,
    notes: 'Vertex and face elements, either byte order; vertex colours and other scalar properties kept as attributes; no faces → point cloud.', load: () => pick(import('./importers/surfaces.js'), 'parsePly') },
  { id: 'off', name: 'Geomview OFF', extensions: ['off'], kinds: ['surface', 'points'], category: 'mesh', status: 'native', via: 'importers/surfaces.js', suites: [],
    notes: 'OFF / COFF / NOFF; polygons triangulated.', load: () => pick(import('./importers/surfaces.js'), 'parseOff') },
  { id: 'gltf', name: 'glTF 2.0 (.gltf)', extensions: ['gltf'], kinds: ['surface', 'polyline', 'points'], category: 'mesh', status: 'native', via: 'importers/gltf.js', suites: [1, 2, 4, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16],
    notes: 'Scene graph flattened with node transforms, y-up rotated to z-up; embedded or side-car buffers. Draco / meshopt compressed primitives are not decoded.', load: () => pick(import('./importers/gltf.js'), 'parseGltf') },
  { id: 'glb', name: 'glTF 2.0 binary (.glb)', extensions: ['glb'], kinds: ['surface', 'polyline', 'points'], category: 'mesh', status: 'native', via: 'importers/gltf.js', suites: [1, 2, 4, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16],
    notes: 'As glTF. Draco / meshopt compressed primitives are not decoded.', load: () => pick(import('./importers/gltf.js'), 'parseGltf') },
  { id: '3mf', name: '3MF', extensions: ['3mf'], kinds: ['surface'], category: 'mesh', status: 'native', via: 'importers/gltf.js', suites: [1, 11, 15, 16],
    notes: 'Core mesh objects, components and build transforms with the declared unit. Beam lattices, slices and materials are not read.', load: () => pick(import('./importers/gltf.js'), 'parse3mf') },

  // ---- tables, centrelines, application format ----
  { id: 'table', name: 'Delimited table (CSV, TSV, XYZ, whitespace)', extensions: ['csv', 'tsv', 'tab', 'xyz', 'txt', 'dat'], kinds: ['points', 'grid', 'polyline', 'frame', 'surface', 'section'], category: 'table', status: 'native', via: 'importers/tables.js', suites: ALL,
    notes: 'Columns mapped to roles by header name (reported in meta.mapping; override with opts.mapping / opts.schema): soundings → points or grid, KP / easting / northing / depth routes → polyline, node + member tables → frame, hull offsets → lofted surface, blade and layout tables, cross-section components.', load: () => pick(import('./importers/tables.js'), 'parseTable') },
  { id: 'route-table', name: 'Route / centreline table (KP, easting, northing, depth)', extensions: ['csv', 'tsv', 'txt', 'xyz', 'dat'], kinds: ['polyline'], category: 'centreline', status: 'native', via: 'importers/tables.js', suites: [8, 9, 10, 12, 13, 14, 15, 16],
    notes: 'The delimited-table reader with the rows taken in order as a line (pipeline, cable, riser, mooring line). One polyline per value of a line / route / branch column; other columns become per-point attributes; attrs.kp is in metres. Chosen with opts.format, or automatically when the headers include KP / chainage.', load: () => pick(import('./importers/tables.js'), 'parseTable', { schema: 'route' }) },
  { id: 'frame-table', name: 'Node and member tables', extensions: ['csv', 'tsv', 'txt', 'dat'], kinds: ['frame', 'points'], category: 'structural', status: 'native', via: 'importers/tables.js', suites: [1, 2, 3, 4, 6, 8, 15],
    notes: 'A node table (id, x, y, z) and a member table (node_a, node_b, D, t or section name, material, group) in one file under [nodes] / [members] headings or as two files; or one table with member end coordinates x1 … z2. Chosen with opts.format, or automatically from the headers.', load: () => pick(import('./importers/tables.js'), 'parseTable', { schema: 'frame' }) },
  { id: 'hull-offsets', name: 'Hull offset table (station, waterline, half-breadth)', extensions: ['csv', 'tsv', 'txt', 'dat'], kinds: ['surface'], category: 'table', status: 'native', via: 'importers/tables.js', suites: [11, 15],
    notes: 'Long form (one offset per row) or matrix form (stations × waterlines): stations lofted linearly and mirrored about the centre plane. Deck, transom and bow are not closed. Chosen with opts.format, or automatically from the headers.', load: () => pick(import('./importers/tables.js'), 'parseTable', { schema: 'hull' }) },
  { id: 'blade-table', name: 'Blade table (radius, chord, twist)', extensions: ['csv', 'tsv', 'txt', 'dat'], kinds: ['polyline'], category: 'table', status: 'native', via: 'importers/tables.js', suites: [16],
    notes: 'Spanwise stations as a polyline along +z with chord (m), twist (deg), relative thickness and airfoil names as attributes. No blade surface is lofted. Chosen with opts.format, or automatically from the headers.', load: () => pick(import('./importers/tables.js'), 'parseTable', { schema: 'blade' }) },
  { id: 'panel-csv', name: 'Panel list (12 coordinates per row)', extensions: ['csv', 'txt', 'dat'], kinds: ['surface'], category: 'panel', status: 'native', via: 'importers/tables.js', suites: [11, 15, 16],
    notes: 'Rows of x1 y1 z1 … x4 y4 z4; quads kept alongside the triangulation.', load: () => pick(import('./importers/tables.js'), 'parseTable', { schema: 'panels' }) },
  { id: 'table-json', name: 'JSON tables (nodes / members / routes / offsets)', extensions: ['json'], kinds: ['points', 'grid', 'polyline', 'frame', 'surface', 'section'], category: 'table', status: 'native', via: 'importers/json.js', suites: [1, 2, 3, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16],
    notes: 'Arrays of records anywhere in the document are read as tables and mapped like CSV columns. There is no standard schema: check meta.mapping.', load: () => pick(import('./importers/json.js'), 'parseTableJson') },
  { id: 'xml-model', name: 'XML tables (nodes / members / routes)', extensions: ['xml'], kinds: ['points', 'grid', 'polyline', 'frame', 'surface', 'section'], category: 'table', status: 'native', via: 'importers/json.js', suites: [1, 2, 3, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16],
    notes: 'Repeated record elements (attributes or simple children) are read as tables and mapped like CSV columns. Schema-specific XML (ISO 15926, vendor formats) is not interpreted beyond that.', load: () => pick(import('./importers/json.js'), 'parseXmlModel') },
  { id: 'oa-json', name: 'Application geometry JSON', extensions: ['json'], kinds: ['surface', 'volume', 'frame', 'polyline', 'grid', 'points', 'section'], category: 'table', status: 'native', via: 'importers/json.js', suites: [],
    notes: 'Lossless round trip of every model kind (written by the JSON exporter).', load: () => pick(import('./importers/json.js'), 'parseOaJson') },

  // ---- finite-element and CFD meshes ----
  { id: 'msh', name: 'Gmsh MSH 2.2 / 4.1 (ASCII)', extensions: ['msh'], kinds: ['frame', 'surface', 'volume', 'points'], category: 'mesh', status: 'native', via: 'importers/gmsh.js', suites: ALL,
    notes: 'First- and second-order elements, physical groups as groups and sets. Binary MSH and MSH 4.0 are rejected with the export setting to use.', load: () => pick(import('./importers/gmsh.js'), 'parseGmsh') },
  { id: 'nastran', name: 'Nastran bulk data', extensions: ['bdf', 'nas', 'nastran'], kinds: ['frame', 'surface', 'volume', 'points'], category: 'structural', status: 'native', via: 'importers/nastran.js', suites: ALL,
    notes: 'GRID, CBAR/CBEAM/CROD/CTUBE/CONROD, CTRIA3/6, CQUAD4/8, CTETRA, CHEXA, CPENTA, CPYRAM; PBAR(L), PBEAM(L), PROD, PTUBE, PSHELL, PSOLID, MAT1; small, large and free field. Loads, constraints and rigid elements are skipped. No unit in the format.', load: () => pick(import('./importers/nastran.js'), 'parseNastran') },
  { id: 'abaqus', name: 'Abaqus input deck', extensions: ['inp'], kinds: ['frame', 'surface', 'volume', 'points'], category: 'structural', status: 'native', via: 'importers/abaqus.js', suites: ALL,
    notes: '*NODE, *ELEMENT (beam, shell, membrane, continuum types), *NSET / *ELSET, beam / shell / solid sections, *MATERIAL, parts and instances. Steps, loads, contact and special elements are skipped. No unit in the format.', load: () => pick(import('./importers/abaqus.js'), 'parseAbaqus') },
  { id: 'ansys-cdb', name: 'ANSYS CDB archive', extensions: ['cdb'], kinds: ['frame', 'surface', 'volume', 'points'], category: 'structural', status: 'native', via: 'importers/ansys.js', suites: [1, 2, 3, 4, 5, 6, 8, 10, 11, 12, 13, 14, 15, 16],
    notes: 'Blocked NBLOCK / EBLOCK (SOLID layout), ET, CMBLOCK, SECTYPE / SECDATA, MP. EBLOCK without the SOLID key and unblocked N / EN commands are not read.', load: () => pick(import('./importers/ansys.js'), 'parseAnsysCdb') },
  { id: 'unv', name: 'I-DEAS universal file', extensions: ['unv'], kinds: ['frame', 'surface', 'volume', 'points'], category: 'mesh', status: 'native', via: 'importers/unv.js', suites: [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16],
    notes: 'Datasets 2411 (nodes), 2412 (elements), 2467 / 2477 (groups), 164 (units). Legacy node datasets 15 / 781 and result datasets are not read.', load: () => pick(import('./importers/unv.js'), 'parseUnv') },
  { id: 'vtk', name: 'VTK legacy (.vtk, ASCII and binary)', extensions: ['vtk'], kinds: ['frame', 'surface', 'volume', 'points', 'grid'], category: 'mesh', status: 'native', via: 'importers/vtk.js', suites: ALL,
    notes: 'Unstructured grid, poly data, structured grid, rectilinear grid, structured points; point-data arrays kept. Cell data and polyhedral cells are not kept.', load: () => pick(import('./importers/vtk.js'), 'parseVtk') },
  { id: 'vtu', name: 'VTK XML (.vtu / .vtp)', extensions: ['vtu', 'vtp'], kinds: ['frame', 'surface', 'volume', 'points'], category: 'mesh', status: 'native', via: 'importers/vtk.js', suites: ALL,
    notes: 'ASCII, inline base64 and appended (base64 / raw) arrays, zlib compression; point-data arrays kept. Image / structured / rectilinear / parallel XML files (.vti, .vts, .vtr, .pvtu) are not read.', load: () => pick(import('./importers/vtk.js'), 'parseVtkXml') },
  { id: 'su2', name: 'SU2 mesh', extensions: ['su2'], kinds: ['frame', 'surface', 'volume'], category: 'mesh', status: 'native', via: 'importers/cfd.js', suites: [7, 8, 11, 15, 16],
    notes: '2-D and 3-D elements with boundary markers as groups (marker edges of a 2-D mesh come back as a frame-kind model named "boundary edges").', load: () => pick(import('./importers/cfd.js'), 'parseSu2') },
  { id: 'openfoam', name: 'OpenFOAM polyMesh', extensions: [], kinds: ['surface'], category: 'mesh', status: 'native', via: 'importers/cfd.js', suites: [7, 8, 11, 15, 16],
    notes: 'ASCII points / faces / owner / boundary, uploaded together or zipped: boundary patches only, as surface groups. The polyhedral cells are not reconstructed; binary polyMesh files are rejected.', load: () => pick(import('./importers/cfd.js'), 'parseOpenFoam') },
  { id: 'exodus', name: 'Exodus II on classic NetCDF', extensions: ['e', 'exo', 'ex2', 'exii', 'g'], kinds: ['frame', 'surface', 'volume', 'points'], category: 'mesh', status: 'native', via: 'importers/netcdf.js', suites: ALL,
    notes: 'Element blocks (bar, tri, quad, tet, hex, wedge, pyramid) and node sets. Side sets and results are not read. Exodus files written on NetCDF-4 / HDF5 need the hdf5 adapter.', load: () => pick(import('./importers/netcdf.js'), 'parseNetcdf') },
  { id: 'exodus-hdf5', name: 'Exodus II on NetCDF-4 / HDF5', extensions: ['e', 'exo', 'ex2', 'exii', 'g'], kinds: ['frame', 'surface', 'volume'], category: 'mesh', status: 'adapter', via: 'adapter:hdf5', suites: ALL,
    notes: 'Read by the HDF5 reader: element blocks (bar, tri, quad, tet, hex, wedge, pyramid) with their names, and node sets. Side sets and results are not read.', missing: 'Convert the file to the classic format, which is read directly: nccopy -k cdf5 in.e out.e', load: adapter('hdf5') },
  { id: 'med', name: 'MED (SALOME)', extensions: ['med', 'rmed'], kinds: ['frame', 'surface', 'volume'], category: 'mesh', status: 'adapter', via: 'adapter:hdf5', suites: ALL,
    notes: 'Read by the HDF5 reader (MED 3 / 4 layout; older files without computation-step groups are accepted): SE2, TR3, QU4, TE4, HE8, PE6, PY5 with the MED node order mapped to this application; quadratic cells on their corner nodes; families and groups become node and element sets. The first mesh and its first computation step are read (opts.mesh chooses another mesh). Polygons, polyhedra, structured MED grids and fields are not read.', missing: 'Export the mesh from SALOME as UNV, or convert it to VTK / Gmsh (for example "meshio convert in.med out.vtu").', load: adapter('hdf5') },
  { id: 'cgns', name: 'CGNS', extensions: ['cgns'], kinds: ['surface', 'volume'], category: 'mesh', status: 'adapter', via: 'adapter:hdf5', suites: [3, 7, 8, 11, 13, 15, 16],
    notes: 'HDF5-based CGNS is read by the HDF5 reader: unstructured sections including MIXED, structured zones as hexahedra (3-D) or quadrilaterals (2-D), sections and boundary conditions as sets; zones are placed side by side without merging interface nodes. NGON_n / NFACE_n polyhedral sections and solution data are not read. ADF-based files must first be converted with "cgnsconvert -h".', missing: 'Convert the mesh to VTK, SU2 or Gmsh (for example "meshio convert in.cgns out.vtu").', load: adapter('hdf5') },
  { id: 'fluent', name: 'Fluent mesh / case (.msh, .cas)', extensions: ['cas'], kinds: ['surface', 'frame'], category: 'mesh', status: 'native', via: 'importers/fluent.js', suites: [7],
    notes: 'Legacy .msh / .cas, ASCII and binary (single or double precision, 32- or 64-bit integers; .gz accepted): nodes, faces (line, triangle, quadrilateral, mixed, polygonal) and zone names. Returns the boundary only — every face with one adjacent cell, wound outwards, as one surface with a group per face zone (name and type); a 2-D mesh gives its boundary edges as a frame-kind model named "boundary edges" with the zones as member groups; a surface mesh without cells is returned whole. Cells are not reconstructed (the file stores them only as faces); interior faces, parent faces of hanging-node meshes, node-less interface faces and every other section (variables, trees, partition, periodic shadows) are counted in meta.skipped. Faces between two cells in non-interior zones (baffles, fans) are left out unless opts.twoSided. Coordinates are taken as metres; a mesh never scaled in Fluent may be in the mesher\'s unit (opts.units). The HDF5 files .msh.h5 / .cas.h5 are read by the HDF5 reader with the same result. Checked on 78 legacy and 9 HDF5 files against independent readings (docs/geometry.md, "Checks on real files"); 64-bit integer sections were met in no real file, and big-endian files are not read. Data files (.dat) hold no mesh.', load: () => pick(import('./importers/fluent.js'), 'parseFluent') },

  // ---- terrain, bathymetry and geospatial ----
  { id: 'esri-ascii', name: 'ESRI ASCII grid', extensions: ['asc', 'agr'], kinds: ['grid'], category: 'geospatial', status: 'native', via: 'importers/grids.js', suites: [3, 4, 5, 7, 8, 9, 10, 12, 13, 15],
    notes: 'Corner or centre registration, NODATA → NaN, optional .prj side-car; block-averaged to opts.maxCells.', load: () => pick(import('./importers/grids.js'), 'parseEsriAscii') },
  { id: 'surfer-grd', name: 'Surfer grid (DSAA / DSBB / DSRR)', extensions: ['grd'], kinds: ['grid'], category: 'geospatial', status: 'native', via: 'importers/grids.js', suites: [],
    notes: 'Surfer 6 text and binary and Surfer 7 binary grids. Other .grd dialects (GMT binary, ESRI binary) are not read.', load: () => pick(import('./importers/grids.js'), 'parseSurferGrd') },
  { id: 'geotiff', name: 'GeoTIFF', extensions: ['tif', 'tiff', 'gtiff'], kinds: ['grid'], category: 'geospatial', status: 'native', via: 'importers/geotiff.js', suites: [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],
    notes: 'Strips and tiles; none / LZW / Deflate / PackBits; predictors; classic and BigTIFF; GeoKeys CRS; GDAL_NODATA; block-averaged to opts.maxCells while decoding. JPEG / LERC / ZSTD compression and rotated rasters are rejected with the gdal command to fix them.', load: () => pick(import('./importers/geotiff.js'), 'parseGeoTiff') },
  { id: 'netcdf', name: 'NetCDF classic (CF grid, GMT grid, UGRID mesh)', extensions: ['nc', 'cdf'], kinds: ['grid', 'points', 'surface'], category: 'geospatial', status: 'native', via: 'importers/netcdf.js', suites: [7, 9, 10, 16],
    notes: 'CDF-1 / CDF-2 / CDF-5. Gridded variable on x/y or lon/lat axes (packing, fill values, positive-down handled); unevenly spaced axes come back as points; UGRID 2-D face meshes as a surface.', load: () => pick(import('./importers/netcdf.js'), 'parseNetcdf') },
  { id: 'netcdf4', name: 'NetCDF-4', extensions: ['nc', 'nc4'], kinds: ['grid', 'surface'], category: 'geospatial', status: 'adapter', via: 'adapter:hdf5', suites: [7, 9, 10, 16],
    notes: 'Read by the HDF5 reader. CF grids (GEBCO elevation, ETOPO z, EMODnet naming; packing, fill values, positive-down handled) are block-averaged to opts.maxCells nodes (default 1.5 million) while being read; unevenly spaced axes come back as points; UGRID 2-D face meshes as a surface. The whole file has to be in memory. Variables compressed with plug-in filters (Blosc, Zstandard) cannot be read.', missing: 'Convert the file to the classic format, which is read directly: nccopy -k classic in.nc out.nc (or ncks -3 in.nc out.nc).', load: adapter('hdf5') },
  { id: 'hdf5', name: 'HDF5 with a documented layout', extensions: ['h5', 'hdf5', 'he5'], kinds: ['grid', 'points', 'surface', 'volume'], category: 'geospatial', status: 'adapter', via: 'adapter:hdf5', suites: [3],
    notes: 'The HDF5 reader recognises NetCDF-4, MED, CGNS and Exodus layouts from the object tree whatever the extension. Other HDF5 files have no geometry schema: nothing is imported and the object tree is listed in meta.tree, unless opts.datasets names the datasets ({ grid, x, y }, { points } or { positions, elements }).', missing: 'Export the datasets to a format read directly: a grid as GeoTIFF, ESRI ASCII or classic NetCDF; points as CSV / XYZ; a mesh as VTK.', load: adapter('hdf5') },
  { id: 'geojson', name: 'GeoJSON', extensions: ['geojson', 'geojsonl', 'ndjson', 'json'], kinds: ['points', 'polyline', 'section'], category: 'geospatial', status: 'native', via: 'importers/vector.js', suites: [7, 8, 9, 10, 12, 13, 14, 15, 16],
    notes: 'Points (with properties), line strings, polygons (flat → section loops, 3-D → closed polylines); longitude/latitude projected to metres; legacy crs member honoured. TopoJSON is not read.', load: () => pick(import('./importers/vector.js'), 'parseGeoJson') },
  { id: 'wkt', name: 'Well-known text (WKT / EWKT)', extensions: ['wkt'], kinds: ['points', 'polyline', 'section', 'surface'], category: 'geospatial', status: 'native', via: 'importers/vector.js', suites: [],
    notes: 'POINT … MULTIPOLYGON, GEOMETRYCOLLECTION, TIN and POLYHEDRALSURFACE; SRID prefix. Curved geometry types are not read.', load: () => pick(import('./importers/vector.js'), 'parseWkt') },
  { id: 'shapefile', name: 'ESRI Shapefile (.shp + .dbf + .prj)', extensions: ['shp'], kinds: ['points', 'polyline', 'section'], category: 'geospatial', status: 'native', via: 'importers/shapefile.js', suites: [7, 8, 13, 15, 16],
    notes: 'Point, MultiPoint, PolyLine, Polygon and their Z / M forms with .dbf attributes and .prj CRS; upload the set together or zipped. MultiPatch records are skipped.', load: () => pick(import('./importers/shapefile.js'), 'parseShapefile') },
  { id: 'landxml', name: 'LandXML', extensions: ['xml'], kinds: ['surface', 'polyline', 'points', 'frame'], category: 'geospatial', status: 'native', via: 'importers/gisxml.js', suites: [3, 4, 5, 7, 8, 10, 13, 15, 16],
    notes: 'TIN surfaces, breaklines, CgPoints, alignments (lines, arcs, spirals approximated, profile PVIs with vertical curves), plan features, pipe networks → frame. Parcels, cross-sections and roadway design are not read.', load: () => pick(import('./importers/gisxml.js'), 'parseLandXml') },
  { id: 'gml', name: 'GML / CityGML / InfraGML geometry', extensions: ['gml', 'citygml', 'xml'], kinds: ['points', 'polyline', 'section', 'surface'], category: 'geospatial', status: 'native', via: 'importers/gisxml.js', suites: [3, 5, 7, 8, 9, 10, 12, 13, 14, 15, 16],
    notes: 'All gml:Point / LineString / Curve / Polygon / surface patches in the document; srsName, srsDimension and latitude-first axis order honoured; 3-D polygons of solids and surfaces triangulated per thematic element. Feature attributes, xlink references and implicit geometry are not read.', load: () => pick(import('./importers/gisxml.js'), 'parseGml') },
  { id: 'kml', name: 'KML / KMZ', extensions: ['kml', 'kmz'], kinds: ['points', 'polyline', 'section'], category: 'geospatial', status: 'native', via: 'importers/gisxml.js', suites: [],
    notes: 'Placemark points, line strings, polygons, multi-geometries and gx:Track, with ExtendedData. Overlays, models and network links are skipped.', load: () => pick(import('./importers/gisxml.js'), 'parseKml') },
  { id: 'gpx', name: 'GPX', extensions: ['gpx'], kinds: ['points', 'polyline'], category: 'geospatial', status: 'native', via: 'importers/gisxml.js', suites: [],
    notes: 'Waypoints, routes and track segments.', load: () => pick(import('./importers/gisxml.js'), 'parseGpx') },
  { id: 'geopackage', name: 'GeoPackage', extensions: ['gpkg', 'sqlite'], kinds: ['points', 'polyline', 'section', 'grid'], category: 'geospatial', status: 'adapter', via: 'adapter:sqlite', suites: [5, 7, 8, 13, 15, 16],
    notes: 'Vector feature tables are read by the built-in SQLite reader (no download): points → points, lines → polylines, flat polygons → sections. Tile pyramids (imagery and tiled gridded coverage) are listed, not decoded; WITHOUT ROWID tables are reported, not read; changes still in a "-wal" side-car file are not seen.', missing: 'Export the layers as GeoJSON or Shapefile (vector) or GeoTIFF (raster), for example "ogr2ogr -f GeoJSON out.geojson in.gpkg layer".', load: adapter('sqlite') },

  // ---- point clouds ----
  { id: 'las', name: 'LAS 1.0–1.4 (uncompressed)', extensions: ['las'], kinds: ['points'], category: 'pointcloud', status: 'native', via: 'importers/pointcloud.js', suites: ALL,
    notes: 'Point formats 0–10; intensity, classification and colour kept; CRS from GeoKey or WKT records; thinned to opts.maxPoints.', load: () => pick(import('./importers/pointcloud.js'), 'parseLas') },
  { id: 'laz', name: 'LAZ (compressed LAS)', extensions: ['laz'], kinds: ['points'], category: 'pointcloud', status: 'adapter', via: 'adapter:laz', suites: ALL,
    notes: 'Decompressed by the LAZ reader, then treated as LAS: point formats 0–3 and 6–8 in chunked LASzip streams, thinned to opts.maxPoints (default 2 million) by keeping every n-th point; intensity, classification and colour kept. Waveform formats (4, 5, 9, 10) and the pre-2011 unchunked stream are not decoded; files of 2 GB or more cannot be opened.', missing: 'Decompress to LAS, which is read directly: laszip -i in.laz -o out.las (or pdal translate in.laz out.las).', load: adapter('laz') },
  { id: 'pts', name: 'Leica PTS', extensions: ['pts'], kinds: ['points'], category: 'pointcloud', status: 'native', via: 'importers/pointcloud.js', suites: [1, 2, 4, 5, 7, 8, 11, 14, 15, 16],
    notes: 'x y z intensity r g b rows with optional point count.', load: () => pick(import('./importers/pointcloud.js'), 'parsePts') },
  { id: 'ptx', name: 'Leica PTX', extensions: ['ptx'], kinds: ['points'], category: 'pointcloud', status: 'native', via: 'importers/pointcloud.js', suites: [],
    notes: 'One or more scans, each placed with its registration matrix; empty returns dropped.', load: () => pick(import('./importers/pointcloud.js'), 'parsePtx') },
  { id: 'xyzrgb', name: 'XYZ with colour / intensity / normals', extensions: ['xyzrgb', 'xyzi', 'xyzn', 'txt', 'xyz', 'asc'], kinds: ['points'], category: 'pointcloud', status: 'native', via: 'importers/pointcloud.js', suites: [1, 2, 4, 5, 7, 8, 11, 14, 15, 16],
    notes: 'Header-less rows of 4, 6, 7, 9 or 10 numbers.', load: () => pick(import('./importers/pointcloud.js'), 'parseXyzRgb') },
  { id: 'e57', name: 'ASTM E57', extensions: ['e57'], kinds: ['points'], category: 'pointcloud', status: 'native', via: 'importers/e57.js', suites: ALL,
    notes: 'ASTM E2807 scans: cartesian or spherical coordinates stored as floats, integers or scaled integers of any bit width, read packet by packet with the page checksums (CRC-32C) verified; points flagged invalid are dropped and counted; the pose of every scan is applied and all scans are merged into one cloud (meta.scans lists name, guid and point count); intensity kept as written, colour scaled to 0–255 by colorLimits; a projected system in coordinateMetadata is reported as the CRS; thinned to opts.maxPoints (default 5 million) by keeping every n-th point. Checked on 69 real files against pye57 (docs/geometry.md, "Checks on real files"). Images, row / column / return indices, time stamps and normals are not read; files that use extension codecs are refused.', load: () => pick(import('./importers/e57.js'), 'parseE57') },

  // ---- drawings, structural and CAD exchange ----
  { id: 'dxf', name: 'AutoCAD DXF (ASCII)', extensions: ['dxf'], kinds: ['frame', 'surface', 'polyline', 'section', 'points'], category: 'cad', status: 'native', via: 'importers/dxf.js', suites: ALL,
    notes: 'LINE → frame, LWPOLYLINE / POLYLINE (bulges, polyface, mesh), 3DFACE, SOLID, CIRCLE, ARC, ELLIPSE, SPLINE, POINT, INSERT; layers; $INSUNITS; closed flat loops → section. 3DSOLID, REGION, BODY and SURFACE entities are read through the ACIS readers (SAT text of R13 – R2010 drawings with its character substitution undone, SAB data of the ACDSDATA section from R2013 on) and give one surface per body — on 158 real drawings (R14 – R2018) 447 of 467 such entities were read, the rest carrying no ACIS data, and all 105 solids came out watertight. Hatches, text and dimensions are not read.', load: () => pick(import('./importers/dxf.js'), 'parseDxf') },
  { id: 'dxf-binary', name: 'AutoCAD binary DXF', extensions: ['dxf'], kinds: ['frame', 'surface', 'polyline', 'section', 'points'], category: 'cad', status: 'native', via: 'importers/dxfbinary.js', suites: [],
    notes: 'The binary group stream (16-bit group codes of R13 and later, single-byte codes with the 255 escape of R12 and earlier) is decoded into the same group pairs as ASCII DXF and read by the same entity handling, so the entities, layers, blocks and units read are exactly those of ASCII DXF. A file that ends before its EOF group is rejected as truncated.', load: () => pick(import('./importers/dxfbinary.js'), 'parseDxfBinary') },
  { id: 'dwg', name: 'AutoCAD DWG', extensions: ['dwg'], kinds: ['frame', 'surface', 'polyline', 'section', 'points'], category: 'cad', status: 'adapter', via: 'adapter:dwg', suites: ALL,
    notes: 'Read with the LibreDWG library (GNU GPL 3), which is NOT included in this application: at the request of the user (opts.consent = true) a separate worker downloads it from the public CDN cdn.jsdelivr.net (one pinned version, about 11 MB) and runs it in the browser; only the drawing data come back, and they are read like a DXF drawing (same entities, layers, blocks and units). Needs a network connection the first time; works in a web browser only. Without consent importFile() throws with the text to show the user. Solids stored as binary ACIS data, hatches, text and dimensions are not read.',
    missing: 'To read DWG, use the application in a web browser with a network connection and confirm the download of the DWG library (LibreDWG, GPL-3.0, not part of this application) when asked. Otherwise save the drawing as ASCII or binary DXF in AutoCAD (or the free ODA File Converter / DWG TrueView); for 3-D solids export STEP.', load: null },
  { id: 'sdnf', name: 'SDNF 3.0', extensions: ['sdnf', 'sdn'], kinds: ['frame', 'surface'], category: 'structural', status: 'native', via: 'importers/sdnf.js', suites: [1, 2, 3, 4, 5, 6, 16],
    notes: 'Packet 10 linear members (section names resolved for pipes, box sections, IPE / HE A / HE B and tabulated W shapes) and Packet 20 plates. Holes, loads, connection details, grids and curved members are not read.', load: () => pick(import('./importers/sdnf.js'), 'parseSdnf') },
  { id: 'cis2', name: 'CIS/2 (CIMsteel) STEP file', extensions: ['stp', 'step', 'p21', 'cis'], kinds: ['frame'], category: 'structural', status: 'native', via: 'importers/step.js', suites: [1, 2, 5, 6],
    notes: 'Analysis-model curve elements on nodes and design-model located prismatic parts → members with resolved section names. Plates, connections, bolts, welds and loads are not read.', load: () => pick(import('./importers/step.js'), 'parseCis2') },
  { id: 'ifc', name: 'IFC 2x3 / 4 / 4.3 — linear elements and alignments', extensions: ['ifc'], kinds: ['frame', 'polyline'], category: 'structural', status: 'native', via: 'importers/ifc.js', suites: [1, 2, 3, 4, 5, 6, 7, 10, 15, 16],
    notes: 'Beams, columns, members, piles and structural curve members → frame with profile, material, GlobalId, name, tag and spatial container; alignment axes → polylines. Surfaces and solids of other products are not produced here (see "IFC — full geometry"). ifcXML and ifcZIP-wrapped ifcXML are not read.', load: () => pick(import('./importers/ifc.js'), 'parseIfc') },
  { id: 'ifc-full', name: 'IFC — full geometry (all products)', extensions: ['ifc'], kinds: ['surface'], category: 'structural', status: 'adapter', via: 'adapter:web-ifc', suites: [1, 2, 3, 4, 5, 6, 7, 10, 15, 16],
    notes: 'Tessellated geometry of every product through the IFC geometry engine (opts.fullGeometry = true selects it): one surface per IFC class with one group per element, openings subtracted, metres, z up. Opening and space volumes only with opts.openings.', missing: 'Without it, beams, columns, members, piles and alignments are still read directly (import without fullGeometry). For plates, slabs and walls export a mesh from the BIM tool (glTF, OBJ or STL).', load: adapter('web-ifc') },
  { id: 'step', name: 'STEP AP203 / AP214 / AP242', extensions: ['step', 'stp', 'p21', 'stpz'], kinds: ['surface', 'polyline', 'points'], category: 'cad', status: 'adapter', via: 'adapter:occt', suites: ALL,
    nativeSubset: 'Tessellated geometry (AP242 triangulated faces and surface sets), faceted B-rep (poly loops), wireframe curves (polylines, lines, circles, B-splines, composite curves) and points are read without the adapter.',
    notes: 'Exact B-rep solids and surfaces need the geometry-kernel reader. Assembly placements are not applied by the built-in subset reader.', missing: 'Without it, export the model as tessellated STEP (AP242 with tessellation), or as a mesh: STL, OBJ or glTF.', load: () => pick(import('./importers/step.js'), 'parseStep') },
  { id: 'iges', name: 'IGES', extensions: ['igs', 'iges'], kinds: ['frame', 'polyline', 'surface', 'points', 'section'], category: 'cad', status: 'adapter', via: 'adapter:occt', suites: ALL,
    nativeSubset: 'Entities 110 (line), 100 (arc), 106 (copious data), 112, 126 (B-spline curve), 102 (composite curve), 116 (point) and independent untrimmed 128 (B-spline surface) are read without the adapter.',
    notes: 'Trimmed and bounded surfaces (144, 143), surfaces of revolution, tabulated cylinders, ruled and offset surfaces and B-rep solids need the geometry-kernel reader.', missing: 'Without it, export wireframe and untrimmed NURBS surfaces (entities 110, 126, 128) or a mesh: STL, OBJ or glTF.', load: () => pick(import('./importers/iges.js'), 'parseIges') },
  { id: 'brep', name: 'OpenCascade BREP', extensions: ['brep', 'brp'], kinds: ['surface'], category: 'cad', status: 'adapter', via: 'adapter:occt', suites: [],
    notes: 'Read and tessellated by the geometry-kernel reader. The format has no unit: millimetres assumed unless stated.', missing: 'Without it, export the shape as tessellated STEP, STL, OBJ or glTF.', load: adapter('occt') },
  { id: 'parasolid', name: 'Parasolid XT (.x_t text, .x_b binary)', extensions: ['x_t', 'x_b', 'xmt_txt', 'xmt_bin', 'xt'], kinds: ['surface', 'frame'], category: 'cad', status: 'native', via: 'importers/parasolid.js', suites: ALL,
    notes: 'Text and binary transmit files (neutral binary, and the bare binary of little-endian machines), the same node stream in either encoding, decoded with the published base schema 13006 and the schema edits newer files embed (Parasolid 14 – 37 as written by SolidWorks, NX, Solid Edge, Onshape and others); files that do not describe their schema (Parasolid 9 – 13, and files saved "as version n") are read with the known layouts of those releases and refused when the stream does not come out consistent — of 125 real files from NX, SolidWorks, IronCAD, Onshape and others (Parasolid 9.0 – 37.0) 112 were read, 146 of their 169 solid bodies came out watertight and 5201 of 5409 faces were read; 10 of the 21 files of Parasolid 9 – 11 were refused. Bodies, multi-body parts and assemblies with instance transforms; faces on PLANE, CYLINDER, CONE, SPHERE, TORUS, B_SURFACE (NURBS) and SWEPT_SURF over a line or B-curve, bounded by LINE, CIRCLE, ELLIPSE, B_CURVE, TRIMMED_CURVE, SP_CURVE and INTERSECTION edges (the latter as the chart of points the file stores), triangulated to opts.deflection (default 1 % of the body diagonal, 15° per facet); one surface per body with one group per face; wire bodies as a frame; coordinates are metres. Binary files: 113 neutral-binary .x_b files written by Solid Edge (Parasolid 36.0) all read, the closed bodies within 0.1 % of the volume of the STEP export of the same part with one exception at 1 % (median 0.02 %, 54 compared); one written by SolidWorks 2009 gives the body of the SolidWorks part it was exported from; one bare-binary file (Parasolid 35.1 saved as 28.1, an assembly of 23 bodies) read with 1 136 of 1 141 faces; and the same decoder reads the partitions inside SolidWorks and Solid Edge part files. Bare binary of big-endian machines is refused. Not read: faces on BLENDED_EDGE (rolling-ball blends kept as a definition — the commonest loss: fillets of SolidWorks and NX parts), OFFSET_SURF, SPUN_SURF and foreign surfaces — counted per type in a warning, with the advice to export STEP; attributes and groups. A node that cannot be interpreted stops the file or the body; nothing is guessed.', load: () => pick(import('./importers/parasolid.js'), 'parseParasolid') },
  { id: 'acis', name: 'ACIS SAT (text)', extensions: ['sat', 'asat'], kinds: ['surface', 'frame'], category: 'cad', status: 'native', via: 'importers/acis.js', suites: ALL,
    notes: 'B-rep of the ACIS / ShapeManager kernel, save versions 1.03 – 28 and the ASM releases of Autodesk programs: faces on plane, cone (cylinders and cones, elliptical ratio), sphere and torus surfaces (apple tori included) and on spline surfaces that store a B-spline (exactsur, the approximation kept with rbblnsur / sumsur / skinsur / cylsur …, or a sum surface built from its two curves), all trimmed by the face boundary; edges on straight, circular and elliptical curves and on intcurves that store their B-spline, cut by the parameters the edge record gives; tolerant topology; "{ ref n }" sub-type references; body transforms and the header unit. Each face is triangulated between its boundary loops to opts.deflection (default: chord error 1 % of the body diagonal, 15° per facet); one surface per body with one group per face; files holding only wire bodies give a frame. Run on 155 real SAT files (Inventor, IronCAD, Revit, Fusion 360, ODA-based writers, ACIS 1.03 – 28): 3789 of 3863 faces read, 135 of 151 solid bodies watertight; volumes of 94 bodies within 0.2 % of the STEP export of the same part (docs/geometry.md, "Checks on real files"). Not read: procedural spline surfaces without a stored B-spline (helical sweeps of threads: sweepsur) — skipped and counted per reason with the advice to export STEP; intcurves without stored data become straight chords (reported); pcurves, attributes and history are passed over.', load: () => pick(import('./importers/acis.js'), 'parseAcis') },
  { id: 'acis-sab', name: 'ACIS SAB (binary)', extensions: ['sab', 'asab'], kinds: ['surface', 'frame'], category: 'cad', status: 'native', via: 'importers/sab.js', suites: [],
    notes: 'The binary form of ACIS SAT: the tagged value stream ("ACIS BinaryFile", "ASM BinaryFile4", and "ASM BinaryFile8" with 64-bit integers) is decoded into the records of the SAT reader, which reads the geometry (same surfaces, curves, options and limits as ACIS SAT). Checked on the SAB data of 159 solids and regions in DXF drawings written by AutoCAD and ODA-based programs (ACIS 218 – 223); the 64-bit variant is implemented from the description of the format and was not met in a real file.', load: () => pick(import('./importers/sab.js'), 'parseSab') },
  { id: 'solidworks', name: 'SolidWorks part (.sldprt)', extensions: ['sldprt', 'sldasm'], kinds: ['surface', 'frame'], category: 'cad', status: 'native', via: 'importers/cadcontainers.js', suites: [2, 9, 11],
    notes: 'The Parasolid bodies of the part are taken out of the file and read by the Parasolid reader, so the geometry read — and not read (rolling-ball blends, offsets …, counted in a warning) — is that of a Parasolid file; metres. Both containers: the compound (OLE2) file written up to SolidWorks 2014 and the record file written since 2015; the bodies are the zlib-compressed Parasolid partition of the stream Contents/Config-n-Partition (lowest n by default, opts.configuration chooses; rollback data and the "ghost" partition are not read). Run on 75 part files from 20 public repositories (64 record files, 11 compound files; Parasolid 20 – 36): 62 read with 105 bodies, 84 of them closed, 15 776 of 16 272 faces; of 17 parts that exist as Parasolid or STEP export of the same name, 13 agree with it in volume within 0.08 % (the seven Parasolid exports of the same revision to six digits), 3 are another revision (different face count) and one multi-body part holds 20 bodies where its STEP export has 12 — bodies hidden in SolidWorks are in the file and are read, with a warning for every multi-body part. Not read: parts made of imported geometry, which keep one unknitted body per import feature (FeatureBodies/LocalBodies) and an empty partition — 10 of the 75 files, refused with that explanation; parts without bodies (3); assemblies (.sldasm), which hold references to part files; features, sketches, configurations without stored bodies, colours and materials.', load: () => pick(import('./importers/cadcontainers.js'), 'parseSolidWorks') },
  { id: 'solid-edge', name: 'Solid Edge part (.par, .psm)', extensions: ['par', 'psm', 'pwd'], kinds: ['surface', 'frame'], category: 'cad', status: 'native', via: 'importers/cadcontainers.js', suites: [],
    notes: 'The Parasolid partition of the part (stream PARASOLID/STREAMn.P_B of the compound file, neutral binary, uncompressed) is read by the Parasolid reader, so the geometry read and not read is that of a Parasolid file; metres. Run on the 113 part files of one public model set written by one Solid Edge release (Parasolid 36.0): all read, and every one gives exactly the faces and the volume of the .x_b file that Solid Edge exported of the same part (113 bodies, 87 closed, 13 229 of 13 654 faces; the closed ones within 0.1 % of the volume of the STEP export, median 0.02 %, one at 1 %). Files of other Solid Edge releases, sheet-metal parts (.psm) and weldments (.pwd) use the same stream but were not available for checking. Assemblies (.asm) hold references to part files and are refused with that explanation; features, sketches, construction bodies kept outside the partition, colours and materials are not read.', load: () => pick(import('./importers/cadcontainers.js'), 'parseSolidEdge') },
  { id: 'rhino-3dm', name: 'Rhino 3DM (openNURBS archive)', extensions: ['3dm'], kinds: ['surface', 'polyline', 'points'], category: 'cad', status: 'native', via: 'importers/rhino3dm.js', suites: [2, 9, 11],
    notes: 'Archive versions 2 – 8 (Rhino 2 to Rhino 8, the openNURBS toolkit and rhino3dm, MoI): chunk tree decoded from the openNURBS description; unit system of the settings table (all openNURBS units, custom units by their scale), layers by full path, object names. ON_Mesh as stored (triangles and quads, double-precision vertices when present); ON_Brep from the render mesh Rhino stores with each face, else the analysis mesh, welded per object; ON_Extrusion from its cached mesh, else from profile and path; ON_LineCurve, ON_ArcCurve, ON_PolylineCurve, ON_NurbsCurve and ON_PolyCurve as polylines; ON_Point per layer and ON_PointCloud as point sets; block instances with their transforms, nested. One surface per layer with one group per object name. A face without a stored mesh (files saved with "Save Small", files of the toolkit) is triangulated here in the parameter plane of its surface (ON_NurbsSurface, ON_PlaneSurface, ON_RevSurface, ON_SumSurface) between its trimming loops, to opts.deflection (default 1 % of the object diagonal, at most 2 % of the smaller extent of the face) and opts.angularStep (15°); by default the facets of one face are bounded in number, and faces that needed more are counted in the warning — their area can then be off by several per cent; give deflection or angularStep to lift the bound. Checked against openNURBS itself (rhino3dm 8.35) on 207 real files, 9919 objects: 2495 meshes and the 6064 face meshes of 1084 B-reps agree in vertex and face counts, bounding box and area, 5440 curves in their end and mid points, units and layers in every file; the triangulation made here agrees with Rhino\'s render meshes in area within 1 % on 13 of 15 files where both exist (4.7 % on one with facets at the bound; one with 2 faces refused). Hidden objects and objects on hidden layers are left out and counted unless opts.hidden; opts.renderMeshes = false ignores the stored meshes. Not read, counted per class in a warning: SubD (convert with ToNURBS or Mesh in Rhino), annotations, text dots, hatches, lights, details, clipping planes, linked blocks, extrusions with mitred ends and no mesh; version 1 archives (Rhino 1) are refused; materials, textures, groups and mesh normals are not kept.', load: () => pick(import('./importers/rhino3dm.js'), 'parseRhino3dm') },
  { id: 'nx', name: 'Siemens NX part (.prt)', extensions: ['prt'], kinds: [], category: 'cad', status: 'convert', via: null, suites: [2, 9, 11],
    notes: 'Not read. The container is known (signature SPLMSSTR, a directory of named sections such as /Root/UG_PART/UG_PART with offsets and sizes), but the bodies are inside the object store of that section in an encoding of NX itself: in the two NX part files examined no Parasolid transmit stream was found either in clear or as a zlib stream, so there is nothing the Parasolid reader could be given. NX exports Parasolid text and binary files, which are read.', convert: 'NX part files are not read. In NX export Parasolid (File → Export → Parasolid, text or binary), which is read directly, or STEP (AP214 / AP242).', load: null },
  { id: 'inventor', name: 'Autodesk Inventor part (.ipt)', extensions: ['ipt', 'iam'], kinds: ['surface', 'frame'], category: 'cad', status: 'native', via: 'importers/cadcontainers.js', suites: [2, 9, 11],
    notes: 'The bodies of the part are a binary ACIS (ShapeManager) stream inside the B-rep segment of the compound file (the stream RSeStorage/B… whose name record says "…BRepSegment", zlib-compressed); it is taken out and read by the ACIS SAB reader, so the geometry read and not read is that of a SAB file; the file is in centimetres, the model in metres. Run on 219 part files from 9 public repositories (Inventor 2013, 2014, 2021, 2024 and 2025; builds 170 – 294): 208 read with 967 bodies (756 of them the sheets of one imported surface model), 191 closed, 11 370 of 11 403 faces; 35 closed parts that exist as SAT export of the same name agree with it in volume to ten digits, one with its STEP export within 0.04 %. Not read: files of Inventor 2026 (builds 300 and later), whose segments are compressed with Zstandard — 8 files, refused with the advice to save as SAT or STEP; parts without bodies (sketches or parameters only, 3 files); assemblies (.iam), presentations and drawings, which hold references; sketches, work geometry and feature history in the other segments; colours and materials.', load: () => pick(import('./importers/cadcontainers.js'), 'parseInventor') },
  { id: 'catia', name: 'CATIA V5 / 3DEXPERIENCE (.CATPart, .CATProduct)', extensions: ['catpart', 'catproduct'], kinds: [], category: 'cad', status: 'convert', via: null, suites: [2, 9, 11],
    notes: 'Not read, and not attempted beyond recognising the file (header V5_CFV2): the geometry is that of the CGM kernel in a proprietary, compressed encoding for which no public description exists.', convert: 'CATIA files are not read. In CATIA save the part as STEP (AP214 / AP242) or IGES; for a mesh, STL.', load: null },
  { id: 'creo', name: 'PTC Creo / Pro/ENGINEER (.prt, .asm)', extensions: ['prt', 'asm'], kinds: [], category: 'cad', status: 'convert', via: null, suites: [2, 9, 11],
    notes: 'Not read, and not attempted beyond recognising the file (text header "#UGC:" followed by binary sections): the geometry is that of the Granite kernel in a proprietary encoding. The extension .prt is shared with NX and .asm with Solid Edge; the content decides.', convert: 'Creo and Pro/ENGINEER files are not read. In Creo use Save a Copy → STEP (AP214 / AP242), or Parasolid (.x_t), both of which are read.', load: null },
  { id: 'jt', name: 'JT (.jt)', extensions: ['jt'], kinds: [], category: 'cad', status: 'convert', via: null, suites: [2, 9, 11],
    notes: 'Not read, and not attempted: JT (ISO 14306) keeps tessellation in LZMA- and arithmetic-coded segments and exact geometry as JT B-rep or embedded Parasolid segments; none of the decoders is implemented here.', convert: 'JT files are not read. Export STEP (AP242) from the originating program or from JT2Go / Teamcenter Visualization; for a mesh, export STL or OBJ.', load: null },
  { id: 'revit', name: 'Autodesk Revit (.rvt, .rfa)', extensions: ['rvt', 'rfa'], kinds: [], category: 'cad', status: 'convert', via: null, suites: [2, 9, 11],
    notes: 'Not read, and not attempted beyond recognising the file (a compound file with the streams BasicFileInfo, Global/…, Partitions/…): the element data are a proprietary serialisation without public description. Revit exports IFC, which is read.', convert: 'Revit files are not read. In Revit export IFC (File → Export → IFC), which is read; for single solids export SAT (ACIS), for drawings DXF.', load: null },
  { id: 'native-cad', name: 'Other native CAD / DCC files (SketchUp, FBX, 3DS, Fusion 360, Navisworks, MicroStation DGN …)', extensions: ['skp', 'fbx', '3ds', 'f3d', 'nwd', 'dgn'], kinds: [], category: 'cad', status: 'convert', via: null, suites: [2, 9, 11],
    notes: 'Not read. SketchUp, FBX, 3DS, Navisworks and DGN were not attempted. Fusion 360 archives (.f3d) were examined: they are ZIP archives whose B-rep blobs (Breps.BlobParts/*.smb and *.smbh) are binary ACIS data that the SAB reader decodes (all 63 blobs of the 10 archives that could be opened; 10 further archives are compressed with Zstandard, ZIP method 93, which is not decoded). But one archive holds a blob per component and per history state, with sketch and construction bodies among them, and which bodies are in the design and where its components are placed is recorded in the undocumented design stream — so nothing is offered rather than a guess.', convert: 'Native CAD files are not read. In the originating program export STEP (AP214 / AP242) for solids, IFC for building and infrastructure models, DXF for drawings, or STL / OBJ / glTF for meshes.', load: null },

  // ---- profiles and hydrodynamic panels ----
  { id: 'svg', name: 'SVG vector profile', extensions: ['svg'], kinds: ['section', 'polyline'], category: 'section', status: 'native', via: 'importers/profiles.js', suites: [12, 13, 14],
    notes: 'Closed paths, rects, circles, ellipses and polygons → section loops (curves flattened, transforms applied, y flipped), grouped by layer. A drawing format: physical size from width/height units, otherwise stated by opts.units.', load: () => pick(import('./importers/profiles.js'), 'parseSvg') },
  { id: 'airfoil', name: 'Airfoil coordinates (Selig / Lednicer .dat)', extensions: ['dat', 'txt', 'air', 'afl'], kinds: ['section'], category: 'airfoil', status: 'native', via: 'importers/profiles.js', suites: [16],
    notes: 'One closed loop in fractions of chord (opts.chord scales it); maximum thickness and camber reported.', load: () => pick(import('./importers/profiles.js'), 'parseAirfoil') },
  { id: 'gdf', name: 'WAMIT GDF (low-order panels)', extensions: ['gdf'], kinds: ['surface'], category: 'panel', status: 'native', via: 'importers/panels.js', suites: [11, 16],
    notes: 'Quadrilateral panels with ISX / ISY symmetry expanded (opts.expandSymmetry = false keeps the stored part); quads kept alongside the triangulation. Higher-order (B-spline patch) GDF is not read.', load: () => pick(import('./importers/panels.js'), 'parseGdf') },
  { id: 'nemoh', name: 'NEMOH mesh', extensions: ['nemoh', 'dat'], kinds: ['surface'], category: 'panel', status: 'native', via: 'importers/panels.js', suites: [11, 16],
    notes: 'Node table and quadrilateral panels with xOz symmetry expanded.', load: () => pick(import('./importers/panels.js'), 'parseNemoh') },
  { id: 'hydrostar', name: 'HydroStar mesh (.hst)', extensions: ['hst'], kinds: ['surface'], category: 'panel', status: 'native', via: 'importers/panels.js', suites: [11],
    notes: 'Node table (COORDINATES … ENDCOORDINATES, with or without node numbers; panels name nodes by their numbers, which need not start at 1 or be consecutive) and panel table (PANEL TYPE 0 and TYPE 1, quadrilaterals and triangles); SYMMETRY / SYMMETRY_BODY 1 (half body, mirrored about y = 0) and 2 (quarter body, mirrored about y = 0 and x = 0) expanded (opts.expandSymmetry = false keeps the stored part); bodies of a multi-body file separated through NUMPANEL, each a group; quads kept alongside the triangulation; metres. RHO, GRAVITY, REFLENGTH, damping zones and other keywords are not geometry and are listed in a warning. The format has no public specification: checked on a hull written by the HydroStar tool chain (2082 panels, quarter body) against an independent reading of the same file, and on meshes written by meshmagick against Capytaine and meshmagick reading them (panel counts equal, areas within 1e-8). The HST exporter writes one body without symmetry.', load: () => pick(import('./importers/panels.js'), 'parseHst') },
]);

const BY_ID = new Map(FORMATS.map((f) => [f.id, f]));

/** Error raised for formats that are not read (status 'convert') or whose adapter is not in the build. */
export class GeometryFormatError extends Error {
  constructor(message, formatId, reason) { super(message); this.name = 'GeometryFormatError'; this.formatId = formatId; this.reason = reason; }
}

const adapterName = (fmt) => fmt.via.slice('adapter:'.length);
// ids the adapters use for the formats they serve, where they differ from the registry ids
const ADAPTER_FORMAT = { 'ifc-full': 'ifc' };

function missingAdapter(fmt, why) {
  const name = adapterName(fmt);
  return new GeometryFormatError(`${fmt.name} needs the "${name}" reader (src/geometry/adapters/${name}.js), which ${why}. ${fmt.missing || ''}`.trim(), fmt.id, 'adapter-missing');
}

/**
 * Resolve the adapter of a format. `opts.adapters[name]` overrides the module (an object with `parse`), and
 * `opts.adapters[name] = null` declares it unavailable. Adapter modules that export `available()` returning
 * false (interface stubs for runtimes that are not vendored) count as not available.
 */
async function loadAdapter(fmt, opts) {
  const name = adapterName(fmt);
  if (opts.adapters && name in opts.adapters) {
    if (!opts.adapters[name]) throw missingAdapter(fmt, 'is not available in this build');
    return opts.adapters[name];
  }
  let mod;
  try { mod = await loaderOf(name)(); } catch (e) {
    const notThere = !e || !e.message || /Cannot find module|Failed to fetch|ERR_MODULE_NOT_FOUND|error loading dynamically imported module|Importing a module script failed/i.test(e.message) || (e && e.code === 'ERR_MODULE_NOT_FOUND');
    throw missingAdapter(fmt, notThere ? 'is not available in this build' : `failed to load (${e.message})`);
  }
  if (!mod || typeof mod.parse !== 'function') throw missingAdapter(fmt, 'was found but does not export parse(input, opts)');
  if (typeof mod.available === 'function' && !mod.available()) throw missingAdapter(fmt, typeof mod.unavailableReason === 'function' ? mod.unavailableReason() : 'is not included in this build');
  return mod;
}

/**
 * Which WebAssembly readers this build really has.
 * @returns {Promise<Object<string, { available: boolean, title: string|null, sizeBytes: number|null, library: object|null }>>}
 */
export async function adapterStatus() {
  const out = {};
  for (const name of Object.keys(ADAPTER_LOADERS)) {
    try {
      const m = await ADAPTER_LOADERS[name]();
      out[name] = { available: typeof m.parse === 'function' && (typeof m.available !== 'function' || !!m.available()), title: m.title ?? null, sizeBytes: m.sizeBytes ?? null, library: m.library ?? null };
    } catch { out[name] = { available: false, title: null, sizeBytes: null, library: null }; }
  }
  return out;
}

const stripGz = (name) => String(name || '').replace(/\.gz$/i, '').replace(/\.stpz$/i, '.stp').replace(/\.ifczip$/i, '.ifc');
const junk = (n) => /(^|\/)(__MACOSX\/|\.DS_Store$|Thumbs\.db$)/.test(n) || n.endsWith('/');

/**
 * Registry entry for a file, from its name and (optionally) its first bytes. Synchronous: archives are
 * identified from their directory listing without being decompressed.
 * @returns {object|null} entry of FORMATS, or null when the file is not recognised
 */
export function formatFor(name, bytes) {
  let d = detectFormat({ name, bytes });
  if (d && d.id === 'gzip') { const id = EXTENSION_MAP[extOf(stripGz(name))]; return id ? BY_ID.get(id) || null : null; }
  if (d && d.id === 'zip') {
    let names = [];
    try { names = zipEntries(bytes).map((e) => e.name).filter((n) => !junk(n)); } catch { return null; }
    const set = detectFileSet(names);
    if (!set) return null;
    if (set.id) return BY_ID.get(set.id) || null;
    d = set.main ? { id: EXTENSION_MAP[extOf(set.main)] } : null;
  }
  return d && d.id ? BY_ID.get(d.id) || null : null;
}

/** Value for the `accept` attribute of a file input: every known extension plus the containers. */
export function acceptString() {
  const s = new Set();
  for (const f of FORMATS) for (const e of f.extensions) s.add(`.${e}`);
  for (const e of ['zip', 'gz', 'prj', 'dbf', 'shx', 'cpg', 'bin', 'ifczip']) s.add(`.${e}`);
  return [...s].sort().join(',');
}

/**
 * Import a geometry file (or file set).
 * @param {{ name: string, text?: string, bytes?: Uint8Array, files?: Map<string, Uint8Array> }} input
 * @param {object} [opts] importer options (units, mapping, maxCells, crs …) plus:
 *   `format` — registry id to use instead of detection; `fullGeometry` — route IFC to the web-ifc adapter;
 *   `adapters` — `{ name: { parse } }` overrides for the WebAssembly readers
 * @returns {Promise<{ format: string, models: object[], warnings: string[], meta: object }>}
 */
export async function importFile(input, opts = {}) {
  if (!input || (input.bytes === undefined && typeof input.text !== 'string' && !(input.files && input.files.size))) throw new Error('importFile: expected { name, bytes } (or text, or a files map).');
  if (opts.format && !BY_ID.has(opts.format)) throw new RangeError(`importFile: unknown format id "${opts.format}" (see FORMATS for the ids).`);
  let cur = { ...input };
  let detected = null;
  const containers = [];
  for (let depth = 0; depth < 4; depth++) {
    // archives are opened first, whether or not the caller named the format of what is inside
    const d = detectFormat(cur);
    if (d && d.id === 'gzip') {
      cur = { ...cur, bytes: await gunzip(cur.bytes), name: stripGz(cur.name) };
      containers.push('gzip');
      continue;
    }
    if (d && d.id === 'zip') {
      const files = await readZip(cur.bytes, { filter: (n) => !junk(n) });
      for (const [k, v] of cur.files || []) if (!files.has(k)) files.set(k, v);
      if (!files.size) throw new Error(`${cur.name || 'archive'}: the ZIP archive is empty.`);
      containers.push('zip');
      cur = { name: cur.name, files };
      continue;
    }
    // a file set: hand the importer its main file, with the companions alongside
    if (d && d.main && cur.files && cur.bytes === undefined && typeof cur.text !== 'string') cur = { name: d.main, bytes: d.id === 'openfoam' ? undefined : cur.files.get(d.main), files: cur.files };
    detected = opts.format ? { id: opts.format, reason: 'chosen by the caller' } : d;
    break;
  }
  if (!detected || !detected.id) {
    const listing = cur.files && cur.bytes === undefined ? ` Files supplied: ${[...cur.files.keys()].slice(0, 12).join(', ')}.` : '';
    throw new GeometryFormatError(`"${input.name || 'file'}" was not recognised as a geometry format this application reads.${listing} See the format matrix for what is supported; rename the file with its proper extension if it is one of those.`, null, 'unrecognised');
  }
  let fmt = BY_ID.get(detected.id);
  if (!fmt) throw new GeometryFormatError(`"${input.name || 'file'}" looks like ${detected.id} (${detected.reason}), which this application does not read.`, detected.id, 'unrecognised');
  if (fmt.id === 'ifc' && opts.fullGeometry) fmt = BY_ID.get('ifc-full');
  if (fmt.status === 'convert') throw new GeometryFormatError(`${cur.name || input.name}: ${fmt.convert}`, fmt.id, 'convert');

  const adapterOpts = { ...opts, format: ADAPTER_FORMAT[fmt.id] || fmt.id };
  let result;
  if (fmt.status === 'adapter' && fmt.nativeSubset) {
    // read what can be read without a kernel; go to the adapter only when the file needs it
    const native = await (await fmt.load()).parse(cur, { ...opts, allowEmpty: true });
    if (!native.meta.needsAdapter) result = native;
    else {
      const what = native.meta.brep && Object.keys(native.meta.brep).length ? Object.entries(native.meta.brep).map(([k, v]) => `${v} ${k}`).join(', ') : native.meta.kernelEntities ? Object.entries(native.meta.kernelEntities).map(([k, v]) => `${v} × ${k}`).join(', ') : '';
      let mod = null, failure = null;
      try { mod = await loadAdapter(fmt, opts); } catch (e) { if (e instanceof GeometryFormatError) failure = e; else throw e; }
      let kernel = null;
      if (mod) {
        try { kernel = await mod.parse(cur, adapterOpts); } catch (e) { if (!native.models.length) throw e; failure = e; }
      }
      if (kernel) {
        // the kernel gives the exact surfaces; curves, points and line models come from the built-in reader
        const extra = native.models.filter((m) => m.kind !== 'surface');
        const base = [];
        for (const k in native.meta) if (!(k in (kernel.meta || {})) && !['counts', 'bbox', 'needsAdapter'].includes(k)) base.push(k);
        const meta = { ...(kernel.meta || {}), via: fmt.via };
        for (const k of base) meta[k] = native.meta[k];
        result = finish(kernel.models.concat(extra), (kernel.warnings || []).concat(native.warnings.filter((w) => !/geometry kernel|geometry-kernel|B-rep face/.test(w))), meta);
      } else if (native.models.length) {
        result = native;
        result.warnings.push(failure instanceof GeometryFormatError
          ? `The "${adapterName(fmt)}" reader is not available in this build, so only the part of the file that can be read without it is shown. ${fmt.missing || ''}`.trim()
          : `The geometry-kernel reader could not process the file (${failure.message}); only the part of the file that can be read without it is shown.`);
      } else throw new GeometryFormatError(`${cur.name || input.name}: the file holds only exact geometry${what ? ` (${what})` : ''} that needs the geometry-kernel reader. ${failure.message}`, fmt.id, 'adapter-missing');
    }
  } else if (fmt.status === 'adapter') {
    const mod = await loadAdapter(fmt, opts);
    result = await mod.parse(cur, adapterOpts);
    result.meta = { ...(result.meta || {}), via: fmt.via };
  } else result = await (await fmt.load()).parse(cur, opts);

  if (!result || !Array.isArray(result.models)) throw new Error(`${fmt.name}: the reader returned no result.`);
  const meta = result.meta || {};
  // an importer may recognise a more specific format than the one it was registered under (Exodus in a NetCDF file)
  const id = ['netcdf', 'table-json', 'gltf', 'glb'].includes(fmt.id) && meta.format && BY_ID.has(meta.format) && BY_ID.get(meta.format).via === fmt.via ? meta.format : fmt.id;
  meta.format = id;
  meta.detected = { id: detected.id, reason: detected.reason, containers };
  meta.fileName = cur.name || input.name || null;
  return { format: id, models: result.models, warnings: result.warnings || [], meta };
}
