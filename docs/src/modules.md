[Table of contents](https://petrkryslucsd.github.io/FinEtools.jl)

# Modules

The FinEtools package consists of many modules which fall into several  categories. The top-level module, `FinEtools`, includes all other modules and exports functions to constitute the public interface. The user is free to generate their own public interface, however. More details are provided [here](rollyourown.html).


- Top-level:
     `FinEtools` is the  top-level module.  For interactive use it is enough to do `using FinEtools`, however in some  cases functions from modules need to be  brought into the scope individually (most importantly, the algorithm modules). This is the ONLY  module that EXPORTS  functions, none of the other modules exports a single function. The entire  public (i. e. exported) interface of the FinEtools package is specified  in the file `FinEtools.jl` (i. e. in the `FinEtools` module). The user is free to specify his or her own set of exported functions from the FinEtools package to create an [ad hoc public interface](rollyourown.html).

- Utilities:
`FTypesModule` (types), `PhysicalUnitModule` (definitions of  numbers with physical units), `AssemblyModule` (assembly of elementwise matrices and vectors),   `CSysModule` (coordinate system module),    `MatrixUtilityModule` (utilities for operations on elementwise matrices), `BoxModule`  (support for working with bounding boxes),  `ForceIntensityModule` (force-intensity module),        `RotationUtilModule` (support for spatial rotations).

- Mesh  entities: 
  `FENodeSetModule`, `FESetModule` (node set and finite element set  types).  

- Mesh Generation: 
   `MeshLineModule`,  `MeshQuadrilateralModule`,   `MeshTriangleModule`,   `MeshTetrahedronModule`,             `TetRemeshingModule`,  `VoxelTetMeshingModule`,     `MeshHexahedronModule`,       `VoxelBoxModule`.  

- Mesh manipulation:  `MeshSelectionModule` (searching of nodes  and elements),  `MeshModificationModule` (mesh boundary, merging  of meshes and nodes, smoothing, partitioning),  `MeshUtilModule` (utilities), `FENodeToFEMapModule` (search structure from nodes to elements).

- Mesh import/export:  `MeshImportModule`,  `MeshExportModule`.

- Fields: 
 `FieldModule`,    `GeneralFieldModule`, `ElementalFieldModule`,    `NodalFieldModule` (modules for representing quantities on the mesh).

- Support for  integration over solids, surfaces, curves, and points: 
 `IntegRuleModule`,   `IntegDomainModule`.

- General algorithms: `AlgoBaseModule` (algorithms), `FEMMBaseModule` (FEM machine for general tasks).

- Heat conduction: `AlgoHeatDiffModule` (algorithms), `FEMMHeatDiffModule`, `FEMMHeatDiffSurfModule`  (FEM machines  to evaluate  the  matrix and vector quantities), `MatHeatDiffModule`  (heat diffusion material)

- Acoustics: `AlgoAcoustModule` (algorithms), `FEMMAcoustModule`, `FEMMAcoustSurfModule` (FEM machines to evaluate the matrix and vector quantities),  `MatAcoustFluidModule` (acoustic fluid material).

- Linear deformation:  `AlgoDeforLinearModule` (algorithms), `DeforModelRedModule`, 
`FEMMDeforLinearBaseModule`,  `FEMMDeforLinearModule`, `FEMMDeforLinearMSModule`,  `FEMMDeforWinklerModule` (FEM machines to evaluate the matrix and vector quantities), `MatDeforModule`, `MatDeforElastIsoModule`, `MatDeforElastOrthoModule` (elastic material models).



