# Real CAD sample files

Files written by CAD systems, used by `tests/geometry/acis.real.test.js`, `sab.test.js`, `dxf.acis.test.js` and
`parasolid.test.js`. Each comes unchanged from a public repository whose licence permits redistribution; the two
`.sab` files are the ACIS data of one solid, copied byte for byte out of the `ACDSDATA` section of the drawing
named (the drawings themselves are larger than a fixture should be).

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

The readers were also run on several hundred further files (see `docs/geometry.md`, "Checks on real files") that
are not stored here because their repositories carry another licence, or none.
