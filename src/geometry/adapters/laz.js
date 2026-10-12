// Compressed LiDAR point clouds (LAZ, LASzip).
//
// LAZ decompression needs the laz-perf runtime, which is not vendored in this build.
// The module keeps the adapter interface so that the importer registry can list the format, explain why it
// cannot be opened and tell the user what to convert the file to. To enable it, vendor the runtime under
// vendor/ as described in docs/vendor.md and replace this file with a full adapter exporting the same names.

export const id = 'laz';
export const title = 'Compressed LiDAR point clouds (LAZ, LASzip)';
export const formats = ['laz'];
export const extensions = ['laz'];
export const library = { name: 'laz-perf', version: null, licence: null };
/** Download size of the runtime [bytes]; 0 because it is not part of this build. */
export const sizeBytes = 0;
/** What to convert a file to so that it can be imported without this adapter. */
export const convertTo = 'uncompressed LAS (las2las -i in.laz -o out.las, or pdal translate in.laz out.las)';

/** True when the runtime is part of this build. */
export function available() { return false; }

/**
 * Import a file of this family. Not available in this build.
 * @param {{ name?: string }} input file
 * @returns {Promise<{ models: object[], warnings: string[], meta: object }>} always rejects
 */
export async function parse(input) {
  throw new Error(`${input && input.name ? `${input.name}: ` : ''}LAZ import is not included in this build — convert to ${convertTo}.`);
}
