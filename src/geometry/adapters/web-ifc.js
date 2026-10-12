// IFC building and structural models (IFC2x3, IFC4).
//
// IFC geometry needs the web-ifc runtime, which is not vendored in this build.
// The module keeps the adapter interface so that the importer registry can list the format, explain why it
// cannot be opened and tell the user what to convert the file to. To enable it, vendor the runtime under
// vendor/ as described in docs/vendor.md and replace this file with a full adapter exporting the same names.

export const id = 'web-ifc';
export const title = 'IFC building and structural models (IFC2x3, IFC4)';
export const formats = ['ifc'];
export const extensions = ['ifc', 'ifczip'];
export const library = { name: 'web-ifc', version: null, licence: null };
/** Download size of the runtime [bytes]; 0 because it is not part of this build. */
export const sizeBytes = 0;
/** What to convert a file to so that it can be imported without this adapter. */
export const convertTo = 'STEP AP214/AP242 (most BIM tools export it), glTF/GLB or OBJ';

/** True when the runtime is part of this build. */
export function available() { return false; }

/**
 * Import a file of this family. Not available in this build.
 * @param {{ name?: string }} input file
 * @returns {Promise<{ models: object[], warnings: string[], meta: object }>} always rejects
 */
export async function parse(input) {
  throw new Error(`${input && input.name ? `${input.name}: ` : ''}IFC import is not included in this build — convert to ${convertTo}.`);
}
