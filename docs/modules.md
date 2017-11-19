[Table of contents](https://petrkryslucsd.github.io/FinEtools.jl)

# Modules

FinEtools modules fall into several  categories. It is usually enough to do `using FinEtools`,
but in some  cases modules need to be  brought into the scope individually.

- Top-level:
     `FinEtools` is the  top-level module.   

- Utilities:
`FTypesModule` (types),    `PhysicalUnitModule` (definitions of  numbers with physical units),       
`AssemblyModule` (assembly of elementwise matrices and vectors),   `CSysModule` (coordinate system module),    `MatrixUtilityModule` (utilities for operations on elementwise matrices), `BoxModule`  (support for working with bounding boxes),  `ForceIntensityModule` (force-intensity module),        `RotationUtilModule` (support for spatial rotations).

- Mesh  entities: 
  `FENodeSetModule`, `FESetModule` (node set and finite element set  types).  

- Mesh Generation: 
   `MeshLineModule`,  `MeshQuadrilateralModule`,   `MeshTriangleModule`,   `MeshTetrahedronModule`,             `TetRemeshingModule`,  `VoxelTetMeshingModule`,     `MeshHexahedronModule`,       `VoxelBoxModule`.  

- Mesh manipulation:  `MeshSelectionModule` (searching of nodes  and elements),  `MeshModificationModule` (mesh boundary, merging  of meshes and nodes, smoothing, partitioning),  `MeshUtilModule` (utilities), `FENodeToFEMapModule` (search structure from nodes to elements).

- Mesh import/export:  `MeshImportModule`,  `MeshExportModule`.

- Fields: 
 `FieldModule`,    `GeneralFieldModule`, `ElementalFieldModule`,    `NodalFieldModule` (modules for representing quantities on the mesh).

- Support for  integration: 
 `IntegRuleModule`,   `IntegDataModule`.

- General algorithms: `AlgoBaseModule` (algorithms), `FEMMBaseModule` (FEM machine for general tasks).

- Heat conduction: `AlgoHeatDiffModule` (algorithms), `FEMMHeatDiffModule`, `FEMMHeatDiffSurfModule`  (FEM machines  to evaluate  the  matrix and vector quantities), `MatHeatDiffModule`  (heat diffusion material)

- Acoustics: `AlgoAcoustModule` (algorithms), `FEMMAcoustModule`, `FEMMAcoustSurfModule` (FEM machines to evaluate the matrix and vector quantities),  `MatAcoustFluidModule` (acoustic fluid material).

- Linear deformation:  `AlgoDeforLinearModule` (algorithms), `DeforModelRedModule`, 
`FEMMDeforLinearBaseModule`,  `FEMMDeforLinearModule`, `FEMMDeforLinearMSModule`,  `FEMMDeforWinklerModule` (FEM machines to evaluate the matrix and vector quantities), `MatDeforModule`, `MatDeforElastIsoModule`, `MatDeforElastOrthoModule` (elastic material models).



