# Functions

## Physical units

```@autodocs
Modules = [FinEtools, FinEtools.PhysicalUnitModule]
Private = true
Order = [:function]
```

## Bounding box functions

```@autodocs
Modules = [FinEtools, FinEtools.BoxModule]
Private = true
Order = [:function]
```

## Coordinate systems

```@autodocs
Modules = [FinEtools, FinEtools.CSysModule]
Private = true
Order = [:function]
```

## Matrix utilities

```@autodocs
Modules = [FinEtools, FinEtools.MatrixUtilityModule]
Private = true
Order = [:function]
```

## Finite element sets

```@autodocs
Modules = [FinEtools, FinEtools.FESetModule]
Private = true
Order = [:function]
```

## Finite element nodes

```@autodocs
Modules = [FinEtools, FinEtools.FENodeSetModule]
Private = true
Order = [:function]
```

## Finite element node-to-element map

```@autodocs
Modules = [FinEtools, FinEtools.FENodeToFEMapModule]
Private = true
Order = [:function]
```

## Selecting nodes and elements

```@autodocs
Modules = [FinEtools, FinEtools.MeshSelectionModule]
Private = true
Order = [:function]
```

## Fields

```@autodocs
Modules = [FinEtools, FinEtools.FieldModule, FinEtools.GeneralFieldModule, FinEtools.NodalFieldModule, FinEtools.ElementalFieldModule]
Private = true
Order = [:function]
```

## Integration rule

```@autodocs
Modules = [FinEtools, FinEtools.IntegRuleModule]
Private = true
Order = [:function]
```

## Integration domain

```@autodocs
Modules = [FinEtools, FinEtools.IntegDomainModule]
Private = true
Order = [:function]
```

## Assembly of matrices and vectors

```@autodocs
Modules = [FinEtools, FinEtools.AssemblyModule]
Private = true
Order = [:function]
```

## Meshing

### Meshing with line elements

```@autodocs
Modules = [FinEtools, FinEtools.MeshLineModule]
Private = true
Order = [:function]
```

### Meshing with triangles

```@autodocs
Modules = [FinEtools, FinEtools.MeshTriangleModule]
Private = true
Order = [:function]
```

### Meshing with quadrilaterals

```@autodocs
Modules = [FinEtools, FinEtools.MeshQuadrilateralModule]
Private = true
Order = [:function]
```

### Meshing with tetrahedra

```@autodocs
Modules = [FinEtools, FinEtools.MeshTetrahedronModule, FinEtools.VoxelBoxModule, FinEtools.VoxelTetMeshingModule]
Private = true
Order = [:function]
```

### Meshing with hexahedra

```@autodocs
Modules = [FinEtools, FinEtools.MeshHexahedronModule]
Private = true
Order = [:function]
```

### Mesh modification

```@autodocs
Modules = [FinEtools, FinEtools.MeshModificationModule]
Private = true
Order = [:function]
```

### Meshing utilities

```@autodocs
Modules = [FinEtools, FinEtools.MeshUtilModule]
Private = true
Order = [:function]
```

### Mesh import/export

```@autodocs
Modules = [FinEtools, FinEtools.MeshImportModule, FinEtools.MeshExportModule]
Private = true
Order = [:function]
```

## Vector-cache utilities

```@autodocs
Modules = [FinEtools, FinEtools.VectorCacheModule]
Private = true
Order = [:function]
```

## Surface-normal utilities

```@autodocs
Modules = [FinEtools, FinEtools.SurfaceNormalModule]
Private = true
Order = [:function]
```

## Force intensity

```@autodocs
Modules = [FinEtools, FinEtools.ForceIntensityModule]
Private = true
Order = [:function]
```

## Rotation utilities

```@autodocs
Modules = [FinEtools, FinEtools.RotationUtilModule]
Private = true
Order = [:function]
```

## FEM machines

### Base

```@autodocs
Modules = [FinEtools, FinEtools.FEMMBaseModule]
Private = true
Order = [:function]
```

### Heat diffusion

```@autodocs
Modules = [FinEtools, FinEtools.FEMMHeatDiffModule, FinEtools.FEMMHeatDiffSurfModule]
Private = true
Order = [:function]
```

### Acoustics

```@autodocs
Modules = [FinEtools, FinEtools.FEMMAcoustModule, FinEtools.FEMMAcoustSurfModule]
Private = true
Order = [:function]
```

### Linear deformation

#### Model reduction types

```@autodocs
Modules = [FinEtools, FinEtools.DeforModelRedModule]
Private = true
Order = [:function]
```

#### Base functionality

```@autodocs
Modules = [FinEtools, FinEtools.FEMMDeforLinearBaseModule, FinEtools.FEMMDeforLinearModule, FinEtools.FEMMDeforWinklerModule, FinEtools.FEMMDeforLinearMSModule, FinEtools.FEMMDeforSurfaceDampingModule, FinEtools.FEMMDeforLinearNICEModule, FinEtools.FEMMDeforLinearESNICEModule]
Private = true
Order = [:function]
```

#### Simple FE models

```@autodocs
Modules = [FinEtools, FinEtools.FEMMDeforLinearModule, FinEtools.FEMMDeforWinklerModule,  FinEtools.FEMMDeforSurfaceDampingModule, ]
Private = true
Order = [:function]
```

#### Advanced FE models

```@autodocs
Modules = [FinEtools, FinEtools.FEMMDeforLinearMSModule, FinEtools.FEMMDeforLinearNICEModule, FinEtools.FEMMDeforLinearESNICEModule]
Private = true
Order = [:function]
```

## Algorithms

### Base

```@autodocs
Modules = [FinEtools, FinEtools.AlgoBaseModule]
Private = true
Order = [:function]
```

### Heat conduction

```@autodocs
Modules = [FinEtools, FinEtools.AlgoHeatDiffModule]
Private = true
Order = [:function]
```

### Acoustics

```@autodocs
Modules = [FinEtools, FinEtools.AlgoAcoustModule]
Private = true
Order = [:function]
```

### Linear deformation

```@autodocs
Modules = [FinEtools, FinEtools.AlgoDeforLinearModule]
Private = true
Order = [:function]
```

## Material models

### Material model abstractions

```@autodocs
Modules = [FinEtools, FinEtools.MatModule]
Private = true
Order = [:function]
```

### Material models for acoustics

```@autodocs
Modules = [FinEtools, FinEtools.MatAcoustFluidModule]
Private = true
Order = [:function]
```

### Material models for heat diffusion

```@autodocs
Modules = [FinEtools, FinEtools.MatHeatDiffModule]
Private = true
Order = [:function]
```

### Material for deformation, base functionality

```@autodocs
Modules = [FinEtools, FinEtools.MatDeforModule]
Private = true
Order = [:function]
```

### Material models for elasticity

```@autodocs
Modules = [FinEtools, FinEtools.MatDeforLinearElasticModule, FinEtools.MatDeforElastIsoModule, FinEtools.MatDeforElastOrthoModule,]
Private = true
Order = [:function]
```
