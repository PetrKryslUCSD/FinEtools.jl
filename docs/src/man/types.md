# Types

## Coordinate systems

```@autodocs
Modules = [FinEtools, FinEtools.CSysModule]
Private = true
Order = [:type]
```

## Finite element sets

```@autodocs
Modules = [FinEtools, FinEtools.FESetModule]
Private = true
Order = [:type]
```

## Finite element nodes

```@autodocs
Modules = [FinEtools, FinEtools.FENodeSetModule]
Private = true
Order = [:type]
```

## Finite element node-to-element map

```@autodocs
Modules = [FinEtools, FinEtools.FENodeToFEMapModule]
Private = true
Order = [:type]
```
 
## Fields

```@autodocs
Modules = [FinEtools, FinEtools.FieldModule, FinEtools.GeneralFieldModule, FinEtools.NodalFieldModule, FinEtools.ElementalFieldModule]
Private = true
Order = [:type]
```

## Integration rule

```@autodocs
Modules = [FinEtools, FinEtools.IntegRuleModule]
Private = true
Order = [:type]
```

## Integration domain

```@autodocs
Modules = [FinEtools, FinEtools.IntegDomainModule]
Private = true
Order = [:type]
```

## Assembly of matrices and vectors

```@autodocs
Modules = [FinEtools, FinEtools.AssemblyModule]
Private = true
Order = [:type]
```

## Mesh import/export

```@autodocs
Modules = [FinEtools, FinEtools.MeshImportModule, FinEtools.MeshExportModule]
Private = true
Order = [:type]
```

## Vector-cache utilities

```@autodocs
Modules = [FinEtools, FinEtools.VectorCacheModule]
Private = true
Order = [:type]
```

## Surface-normal utilities

```@autodocs
Modules = [FinEtools, FinEtools.SurfaceNormalModule]
Private = true
Order = [:type]
```

## Force intensity

```@autodocs
Modules = [FinEtools, FinEtools.ForceIntensityModule]
Private = true
Order = [:type]
```

## FEM machines
### Base

```@autodocs
Modules = [FinEtools, FinEtools.FEMMBaseModule]
Private = true
Order = [:type]
```

### Heat diffusion

```@autodocs
Modules = [FinEtools, FinEtools.FEMMHeatDiffModule, FinEtools.FEMMHeatDiffSurfModule]
Private = true
Order = [:type]
```

### Acoustics

```@autodocs
Modules = [FinEtools, FinEtools.FEMMAcoustModule, FinEtools.FEMMAcoustSurfModule]
Private = true
Order = [:type]
```

### Linear deformation

```@autodocs
Modules = [FinEtools, FinEtools.DeforModelRedModule, FinEtools.FEMMDeforLinearBaseModule, FinEtools.FEMMDeforLinearModule, FinEtools.FEMMDeforWinklerModule, FinEtools.FEMMDeforLinearMSModule, FinEtools.FEMMDeforSurfaceDampingModule, FinEtools.FEMMDeforLinearNICEModule, FinEtools.FEMMDeforLinearESNICEModule]
Private = true
Order = [:type]
```

## Material models

### Material model abstractions

```@autodocs
Modules = [FinEtools, FinEtools.MatModule]
Private = true
Order = [:type]
```

### Material models for acoustics

```@autodocs
Modules = [FinEtools, FinEtools.MatModule, FinEtools.MatAcoustFluidModule]
Private = true
Order = [:type]
```

### Material models for heat diffusion

```@autodocs
Modules = [FinEtools, FinEtools.MatModule, FinEtools.MatHeatDiffModule]
Private = true
Order = [:type]
```

### Material for deformation, base functionality

```@autodocs
Modules = [FinEtools, FinEtools.MatDeforModule]
Private = true
Order = [:type]
```

### Material models for elasticity

```@autodocs
Modules = [FinEtools, FinEtools.MatDeforLinearElasticModule, FinEtools.MatDeforElastIsoModule, FinEtools.MatDeforElastOrthoModule,]
Private = true
Order = [:type]
```
