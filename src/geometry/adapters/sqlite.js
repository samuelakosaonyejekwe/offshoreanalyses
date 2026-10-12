// SQLite containers: OGC GeoPackage vector features and tiled gridded coverage.
//
// GeoPackage needs the sql.js runtime, which is not vendored in this build.
// The module keeps the adapter interface so that the importer registry can list the format, explain why it
// cannot be opened and tell the user what to convert the file to. To enable it, vendor the runtime under
// vendor/ as described in docs/vendor.md and replace this file with a full adapter exporting the same names.

export const id = 'sqlite';
export const title = 'SQLite containers: OGC GeoPackage vector features and tiled gridded coverage';
export const formats = ['geopackage'];
export const extensions = ['gpkg', 'sqlite'];
export const library = { name: 'sql.js', version: null, licence: null };
/** Download size of the runtime [bytes]; 0 because it is not part of this build. */
export const sizeBytes = 0;
/** What to convert a file to so that it can be imported without this adapter. */
export const convertTo = 'GeoJSON or Shapefile for vector layers (ogr2ogr -f GeoJSON out.geojson in.gpkg), GeoTIFF for rasters (gdal_translate in.gpkg out.tif)';

/** True when the runtime is part of this build. */
export function available() { return false; }

/**
 * Import a file of this family. Not available in this build.
 * @param {{ name?: string }} input file
 * @returns {Promise<{ models: object[], warnings: string[], meta: object }>} always rejects
 */
export async function parse(input) {
  throw new Error(`${input && input.name ? `${input.name}: ` : ''}GeoPackage / SQLite import is not included in this build — convert to ${convertTo}.`);
}
