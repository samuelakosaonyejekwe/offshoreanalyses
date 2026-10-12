// HDF5 containers: NetCDF-4 grids (CF conventions), MED (Salome / Code_Aster) and CGNS meshes.
//
// HDF5 needs the h5wasm runtime, which is not vendored in this build.
// The module keeps the adapter interface so that the importer registry can list the format, explain why it
// cannot be opened and tell the user what to convert the file to. To enable it, vendor the runtime under
// vendor/ as described in docs/vendor.md and replace this file with a full adapter exporting the same names.

export const id = 'hdf5';
export const title = 'HDF5 containers: NetCDF-4 grids (CF conventions), MED (Salome / Code_Aster) and CGNS meshes';
export const formats = ['hdf5', 'netcdf4', 'med', 'cgns'];
export const extensions = ['h5', 'hdf5', 'he5', 'nc', 'nc4', 'med', 'rmed', 'cgns'];
export const library = { name: 'h5wasm', version: null, licence: null };
/** Download size of the runtime [bytes]; 0 because it is not part of this build. */
export const sizeBytes = 0;
/** What to convert a file to so that it can be imported without this adapter. */
export const convertTo = 'NetCDF classic (nccopy -k classic in.nc out.nc), or for meshes Gmsh .msh, VTK .vtu or CGNS-ADF/Nastran/Abaqus decks';

/** True when the runtime is part of this build. */
export function available() { return false; }

/**
 * Import a file of this family. Not available in this build.
 * @param {{ name?: string }} input file
 * @returns {Promise<{ models: object[], warnings: string[], meta: object }>} always rejects
 */
export async function parse(input) {
  throw new Error(`${input && input.name ? `${input.name}: ` : ''}HDF5 (NetCDF-4, MED, CGNS) import is not included in this build — convert to ${convertTo}.`);
}
