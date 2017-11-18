[Table of contents](https://petrkryslucsd.github.io/FinEtools.jl)

# Modules

FinEtools modules fall into several  categories. It is usually enough to do `using FinEtools`,
but in some  cases modules need to be  brought into the scope individually.

- Top-level:
     `FinEtools` is the  top-level module.   

- Algorithms: The various algorithms are implemented in modules for common methods
`AlgoBaseModule`, for acoustics `AlgoAcoustModule`, for linear deformation `AlgoDeforLinearModule`, and for heat conduction `AlgoHeatDiffModule`.

- Geometry  Data: 
`CSysModule`, `IntegDataModule`.

- Finite Element  Model Machines
   + Base
   `FEMMBaseModule`,
   + Heat diffusion
    `FEMMHeatDiffModule`, `FEMMHeatDiffSurfModule`  
   + Acoustics
    `FEMMAcoustModule`, `FEMMAcoustSurfModule`
   + Linear deformation
    `FEMMDeforLinearBaseModule`, `DeforModelRedModule`,  `FEMMDeforLinearModule`, `FEMMDeforLinearMSModule`,  `FEMMDeforWinklerModule`    

- Materials
  `MatDeforModule`     `MatDeforElastIsoModule`  `MatDeforElastOrthoModule`  `MatHeatDiffModule`  
      `MatAcoustFluidModule`         
  
- Mesh  entities
  `FENodeSetModule`, `FESetModule`   

- Mesh Generation and Manipulation
   `MeshLineModule`,  `MeshQuadrilateralModule`  `MeshTriangleModule`  `MeshTetrahedronModule` `MeshImportModule`     `MeshSelectionModule`         `MeshModificationModule`   `VoxelTetMeshingModule` `MeshExportModule`    `MeshHexahedronModule`  `TetRemeshingModule`   `MeshUtilModule` `FENodeToFEMapModule` 
- Fields 
 `FieldModule`    `GeneralFieldModule` `ElementalFieldModule`    `NodalFieldModule`

- Utilities
`FTypesModule`,    `PhysicalUnitModule`       
`AssemblyModule`      `MatrixUtilityModule` `BoxModule`   `IntegRuleModule`    
 `ForceIntensityModule`        `VoxelBoxModule`    `RotationUtilModule`
         
      
                    
            
            
           
           

