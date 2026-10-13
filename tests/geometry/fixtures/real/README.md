# Real sample files

Files written by CAD systems, scanners' software and meshers, used by `tests/geometry/acis.real.test.js`,
`sab.test.js`, `dxf.acis.test.js`, `parasolid.test.js`, `parasolid.binary.test.js`, `cadcontainers.test.js`,
`e57.test.js`, `fluent.test.js`, `rhino3dm.test.js` and `hydrostar.test.js`. Each comes unchanged from a public
repository whose licence permits redistribution, with three kinds of exception that are named in the table: the two
`.sab` files are the ACIS data of one solid, copied byte for byte out of the `ACDSDATA` section of the drawing named
(the drawings themselves are larger than a fixture should be); the two `inventor-pin.*.bin` files are two streams of the
part file named, copied byte for byte out of its compound-file container for the same reason; and
`meshmagick-box.hst` was written here with meshmagick, since no HydroStar mesh was found under a licence that
allows redistribution.

| File | Bytes | Written by | Source (repository @ commit, path) | Licence of the repository |
| --- | ---: | --- | --- | --- |
| `fusion-box.sat` | 4728 | Autodesk Fusion 360 (ASM 222, SAT 7.0) | <https://github.com/opendesk/opendesk-on-demand> @ `73e6bbf91b2e92bb8ef974dc2e07a1049ccadeec`, `docs/fusion_files/box v1.sat` | Unlicense |
| `acis28-rod.sat` | 3256 | Ansys Electronics Desktop (ACIS 28.0, SAT 28.0) | <https://github.com/ansys/example-data> @ `8ef505f6ef699eeec68a16001df91ba353342f40`, `pyaedt/mri/implant_rod.sat` | MIT |
| `revit-solids.sat` | 77201 | Autodesk Revit 2014 (ASM 219, SAT 7.0) | <https://github.com/jeremytammik/tbc> @ `7a00cb7fbb97690c5e5d94f91f543f101c9c136d`, `a/zip/SolidFile.sat` | MIT |
| `ezdxf-3dsolids-r2010.dxf` | 161075 | ODA-based writer ("Open Design Alliance ACIS Builder"; DXF R2010, 3DSOLID with SAT 7.0 text) | <https://github.com/mozman/ezdxf> @ `7c64bce52c9d60b9db2dad4d2793c1168c8d3f64`, `examples/acistools/3dsolids.dxf` | MIT |
| `ezdxf-sphere-r2018.dxf` | 71793 | ODA-based writer (DXF R2018, 3DSOLID with SAB data) | <https://github.com/mozman/ezdxf> @ `7c64bce52c9d60b9db2dad4d2793c1168c8d3f64`, `examples/xref/sphere.dxf` | MIT |
| `autocad-flange.sab` | 14727 | Autodesk AutoCAD (SAB, save version 218) | <https://github.com/Adrian-Amoddio/wildlife-detection-camera-system> @ `4e3b064b87c3ac75f2380077e099e57d5ff0caa2`, `hardware/pmma-flange-files/PMMA_Flange_v1.dxf` (ACDSRECORD of the 3DSOLID with handle 84) | MIT |
| `autocad-collar.sab` | 7795 | Autodesk AutoCAD (SAB, save version 218) | <https://github.com/abhi-243/SAR-Rover> @ `331476c85b89e633603a1d05a0849f0a242a7e07`, `CAD/collar.dxf` (ACDSRECORD of the 3DSOLID with handle 84) | MIT |
| `onshape-sphere-in-box.x_t` | 10294 | Onshape (Parasolid 33.0, schema 32001) | <https://github.com/SimScaleGmbH/simscale-python-sdk> @ `2989f25489ed77db3e8e1778c1dd24f4a84a4303`, `examples/fixtures/flow-around-sphere-v3.x_t` | MIT |
| `solidworks-bars-v10.x_t` | 17340 | SolidWorks (Parasolid 10.0, schema 10004) | <https://github.com/ansys/example-data> @ `8ef505f6ef699eeec68a16001df91ba353342f40`, `pymechanical/embedding/LONGBAR.x_t` | MIT |
| `solidworks-load-cell.x_t` | 38246 | SolidWorks (Parasolid 21.0, schema 20000) | <https://github.com/psas/sw-cad-common> @ `8cfa9724f87e0628cea148a30ed02fb16915c307`, `Misc/Pull_Tester_Crane/Load_Cell_s.x_t` | BSD-2-Clause |
| `solidworks-plate.x_t` | 19830 | SolidWorks (Parasolid 19.1, schema 19008) | <https://github.com/psas/sw-cad-common> @ `8cfa9724f87e0628cea148a30ed02fb16915c307`, `Misc/Custom_Tools/Grinder_Mount/Horozontal_Flat_Plate.x_t` | BSD-2-Clause |
| `solidworks-bell-crank.x_b` | 20366 | SolidWorks 2009 (Parasolid 19.1 neutral binary, schema 19008) | <https://github.com/psas/sw-cad-airframe-rollcontrol> @ `2e18e7c42ec5895adc62a32ab6506f3d94b6d899`, `Bell_Crank.x_b` | BSD-2-Clause |
| `solidworks-ole-washer.sldprt` | 83968 | SolidWorks (compound file as written up to 2014; Parasolid 24.0 partition, schema 20000) | <https://github.com/cartheur/aiventure-david> @ `bf8a9cdf3869ac4ee951091731026df5a95fbc22`, `design/motors/bearing/AS1226.SLDPRT` | MIT |
| `solidworks-record-sleeve.sldprt` | 89481 | SolidWorks (record file as written from 2015 on; Parasolid 34.0 partition, schema 34101) | <https://github.com/gearfinity/gearfinity> @ `2260c4d942eb7b3cd6f0c6983743c6e998cd9499`, `_all_parts/planet_gear_sleeve_bearing-1.SLDPRT` | MIT |
| `inventor-pin.B.bin`, `inventor-pin.M.bin` | 2865, 991 | Autodesk Inventor (build 180, the 2014 release; streams `RSeStorage/Bwjduekaj5efbff5eenfnwjbxfb` — the B-rep segment, zlib — and `RSeStorage/Mwjduekaj5efbff5eenfnwjbxfb` — its name record) | <https://github.com/bad-alloc-heavy-industries/MXKeyboard> @ `a9116bf41f271594d8e464b5f3577af30fcbcde6`, `case/cirpin.ipt` | BSD-3-Clause |
| `meshmagick-box.hst` | 11046 | meshmagick 3.4 (`mmio.write_HST`) from a Capytaine 3.0 mesh (`mesh_parallelepiped(size=(4, 2, 1), center=(0, 0, -0.5), resolution=(6, 4, 3))`) | generated for this repository; no third-party content | — |
| `e57-pdal-two-scans.e57` | 6144 | CloudCompare (two scans; float coordinates, intensity, 8-bit colour) | <https://github.com/PDAL/PDAL> @ `64c11506cccdf790b9677458d5ee884bdf4e754b`, `test/data/e57/A_moved_B.e57` | BSD-3-Clause |
| `e57-las2e57-colour.e57` | 9216 | LAStools las2e57 (scaled-integer coordinates, 16-bit colour, LAS extension fields) | <https://github.com/asmaloney/libE57Format-test-data> @ `24201531715945ba7d49af6dfcb8e07b6b9b6b7f`, `3rdParty/las2e57/ColourRepresentation.e57` | CC0-1.0 |
| `e57-bad-crc.e57` | 1024 | libE57Format (a page with a wrong checksum) | <https://github.com/asmaloney/libE57Format-test-data> @ `24201531715945ba7d49af6dfcb8e07b6b9b6b7f`, `self/bad-crc.e57` | CC0-1.0 |
| `e57-spherical.e57` | 10240 | libE57Format 3.0.0 (spherical coordinates, half the points flagged invalid) | <https://github.com/dancergraham/e57> (fork of cry-inc/e57) @ `7a7498f679b30588dc9298beb7aafab2245a2d0c`, `testdata/tiny_spherical.e57` | MIT |
| `e57-scaled-intensity.e57` | 2048 | libE57Format 3.0.0 (ScaledInteger intensity and limits) | same repository @ same commit, `testdata/scaled_integer_intensity.e57` | MIT |
| `e57-images.e57` | 16384 | libE57Format 3.0.0 (two points, four images2D) | same repository @ same commit, `testdata/tiny_pc_and_images.e57` | MIT |
| `e57-cloudcompare-rgb.e57` | 34816 | CloudCompare (2090 points subsampled from a Faro `.fls` scan; float coordinates, colour) | same repository @ same commit, `testdata/tinyCartesianFloatRgb.e57` | MIT |
| `e57-pye57-translation.e57` | 3072 | pye57 0.4.17 (pose with a translation and no rotation) | <https://github.com/davidcaron/pye57> @ `b8d20999f3037c580bf1fe1dcb239c50a8f764ad`, `tests/test_data/translationWithoutRotation.e57` | MIT |
| `e57-normals.e57` | 3072 | libE57Format 3.1.1 (surface-normals extension fields) | <https://github.com/davidcaron/pye57> @ `b8d20999f3037c580bf1fe1dcb239c50a8f764ad`, `tests/test_data/testWithNormals.e57` | MIT |
| `fluent-cube-5tets-ascii.cas` | 635 | PETSc test mesh written by hand (header `cube_5tets test mesh`; ASCII, mixed face rows) | <https://github.com/petsc/petsc> @ `3a51d72243acc60ae6dc87e3ff055a255cdf3d78`, `share/petsc/datafiles/meshes/cube_5tets_ascii.cas` | BSD-2-Clause |
| `fluent-cube-5tets-binary.cas` | 833 | the same mesh in binary sections 2010 / 2012 / 2013 (single precision, 32-bit integers) | same repository and commit, `share/petsc/datafiles/meshes/cube_5tets.cas` | BSD-2-Clause |
| `fluent-cubit-square.cas` | 6464 | CUBIT ("CUBIT to Fluent File"; 2-D triangles, ASCII) | same repository and commit, `share/petsc/datafiles/meshes/square.cas` | BSD-2-Clause |
| `fluent-wedge-cylinder.cas` | 13658 | Ansys Fluent 24.2.0 (ASCII case file; wedges and hexahedra) | same repository and commit, `share/petsc/datafiles/meshes/wedge_cylinder.cas` | BSD-2-Clause |
| `fluent-hex-binary.msh` | 26946 | Ansys "Fluent Interface 23.1" (binary double precision, 216 hexahedra) | <https://github.com/ansys/example-data> @ `8ef505f6ef699eeec68a16001df91ba353342f40`, `pyprimemesh/unit_test_examples/hex.msh` | MIT |
| `fluent-species-2d.cas` | 41213 | Ansys "Fluent Interface 23.1" (2-D, binary double precision, four boundary zones) | same repository and commit, `pyfluent/2d_box/pro_species.cas` | MIT |
| `rhino5-mesh-brep-polylines.3dm` | 109146 | Rhinoceros 5.0 (openNURBS 201705225, archive version 50) | <https://github.com/BlockResearchGroup/compas_tna> @ `df01b179222fc297225c35c709cbef900ddedff4`, `data/tutorial/rhinogeometry.3dm` | MIT |
| `rhino3dm-revolved-sphere.3dm` | 8695 | rhino3dm (.NET binding of openNURBS, toolkit version 2415953804, macOS; archive version 80) | <https://github.com/mcneel/rhino-developer-samples> @ `4863b4c4ba281e07b539bd812e8660c092268af6`, `rhino3dm/cs/SampleCSPBRMaterial/fileWithMaterial.3dm` | the MIT permission notice of Robert McNeel & Associates (`License.md` of the repository) |

The readers were also run on several hundred further files (see `docs/geometry.md`, "Checks on real files") that
are not stored here because their repositories carry another licence, or none.
