"""
FinEtools (C) 2017, Petr Krysl

Finite Element tools.  Julia implementation  of the finite element method
for continuum mechanics.
"""
module FinEtools

include("FTypesModule.jl")
using FinEtools.FTypesModule
export FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec

include("MatrixUtilityModule.jl")

include("PhysicalUnitModule.jl")
using FinEtools.PhysicalUnitModule
phun = FinEtools.PhysicalUnitModule.phun
export phun

include("RotationUtilModule.jl")
using FinEtools.RotationUtilModule
export rotmat, skewmat!, rotmat3!

include("CSysModule.jl")
using FinEtools.CSysModule
export CSys
export updatecsmat!

include("FESetModule.jl")
using FinEtools.FESetModule
export FESet, FESet0Manifold, FESet1Manifold, FESet2Manifold, FESet3Manifold
export manifdim
export nodesperelem
export count
export getconn!
export setotherdimension!
export setlabel!
export subset
export cat
export updateconn!
export FESetP1
export FESetL2
export FESetL3
export FESetT3
export FESetQ4
export FESetQ9
export FESetQ8
export FESetT6
export FESetH8
export FESetH20
export FESetH27
export FESetT4
export FESetT10


include("FENodeSetModule.jl")
using FinEtools.FENodeSetModule
export FENodeSet
export spacedim
export xyz3
export count

include("FENodeToFEMapModule.jl")
using FinEtools.FENodeToFEMapModule
export FENodeToFEMap

include("FieldModule.jl")
using FinEtools.FieldModule
export Field
export ndofs
export nents
export gathersysvec
export gathersysvec!
export gathervalues_asvec!
export gathervalues_asmat!
export gatherdofnums!
export numberdofs!
export setebc!
export applyebc!
export scattersysvec!

include("GeneralFieldModule.jl")
using FinEtools.GeneralFieldModule
export GeneralField

include("NodalFieldModule.jl")
using FinEtools.NodalFieldModule
export NodalField
export nnodes

include("ElementalFieldModule.jl")
using FinEtools.ElementalFieldModule
export ElementalField
export nelems

include("MeshUtilModule.jl")
using FinEtools.MeshUtilModule

include("MeshSelectionModule.jl")
using FinEtools.MeshSelectionModule
export connectedelems
export connectednodes
export selectnode
export selectelem
export findunconnnodes

include("MeshExportModule.jl")
using FinEtools.MeshExportModule
export vtkexportmesh
export finealemesh
export graphcsv

include("MeshImportModule.jl")
using FinEtools.MeshImportModule

include("MeshModificationModule.jl")
using FinEtools.MeshModificationModule
export meshboundary
export fusenodes
export compactnodes
export mergemeshes
export mergenmeshes!
export mergenodes
export renumberconn!
export meshsmoothing
export mirrormesh

include("MeshQuadrilateralModule.jl")
using FinEtools.MeshQuadrilateralModule
export Q4annulus
export Q8annulus
export Q4quadrilateral
export Q4elliphole
export Q4block
export Q4blockx
export Q8block
export Q4toQ8
export Q4refine

include("MeshLineModule.jl")
using FinEtools.MeshLineModule
export L2block
export L2blockx

include("MeshTriangleModule.jl")
using FinEtools.MeshTriangleModule
export T3blockx
export T3blockx
export T3block
export T3toT6
export T6block
export Q4toT3
export T3refine
export Q4toT3

include("MeshHexahedronModule.jl")
using FinEtools.MeshHexahedronModule
export H8block
export H8blockx
export H8sphere
export H8refine
export H8toH27
export H8hexahedron
export H27block
export H8extrudeQ4
export H8spheren
export H20block
export H8toH20
export H8voximg
export H8compositeplatex

include("MeshTetrahedronModule.jl")
using FinEtools.MeshTetrahedronModule
export T4blocka
export T4blockb
export T4blockca
export T4blockcb
export T4block
export T4blockx
export T4toT10
export T10block

include("ForceIntensityModule.jl")
using FinEtools.ForceIntensityModule
export ForceIntensity
export getforce!

include("MatHeatDiffModule.jl")
using FinEtools.MatHeatDiffModule
export MatHeatDiff

include("MatAcoustFluidModule.jl")
using FinEtools.MatAcoustFluidModule
export MatAcoustFluid

include("DeforModelRedModule.jl")
using FinEtools.DeforModelRedModule
export DeforModelRed
export DeforModelRed1D
export DeforModelRed2DStrain
export DeforModelRed2DStress
export DeforModelRed2DAxisymm
export DeforModelRed3D
export nstrains
export stresscomponentmap
export Blmat!

include("MatDeforModule.jl")
using FinEtools.MatDeforModule
export MatDefor
export strain2x2tto3v!
export strain3vto2x2t!
export strain3x3tto6v!
export strain6vto3x3t!
export stress2x2to3v!
export stress3vto2x2t!
export stress3vto3x3t!
export stress4vto3x3t!
export stress6vto3x3t!
export stress3x3tto6v!
export strain9vto6v!
export strain6vto9v!
export stress9vto6v!
export rotstressvec
export strainvectorrotation
export rotatestiffness!
export rotatecompliance!

include("MatDeforElastIsoModule.jl")
using FinEtools.MatDeforElastIsoModule
export MatDeforElastIso

include("MatDeforElastOrthoModule.jl")
using FinEtools.MatDeforElastOrthoModule
export MatDeforElastOrtho

include("AssemblyModule.jl")
using FinEtools.AssemblyModule
export SysmatAssemblerBase
export SysmatAssemblerSparse
export startassembly!
export assemble!
export makematrix!
export SysmatAssemblerSparseSymm
export startassembly!
export assemble!
export makematrix!
export SysvecAssemblerBase
export SysvecAssembler
export startassembly!
export assemble!
export makevector!

include("IntegRuleModule.jl")
using FinEtools.IntegRuleModule
export IntegRule
export TriRule
export GaussRule
export TetRule

include("GeoDModule.jl")
using FinEtools.GeoDModule
export GeoD
export otherdimensionunity
export Jacobianpoint
export Jacobiancurve
export Jacobiansurface
export Jacobianvolume
export Jacobianmdim

include("FEMMBaseModule.jl")
using FinEtools.FEMMBaseModule
export FEMMAbstractBase
export FEMMBase
export integrationdata
export integratefieldfunction
export integratefunction
export distribloads
export fieldfromintegpoints
export elemfieldfromintegpoints

include("FEMMHeatDiffModule.jl")
using FinEtools.FEMMHeatDiffModule
export FEMMHeatDiff
export conductivity
export nzebcloadsconductivity

include("FEMMHeatDiffSurfModule.jl")
using FinEtools.FEMMHeatDiffSurfModule
export FEMMHeatDiffSurf
export surfacetransfer
export surfacetransferloads
export nzebcsurfacetransferloads

include("FEMMAcoustModule.jl")
using FinEtools.FEMMAcoustModule
export FEMMAcoust
export acousticmass
export nzebcloadsacousticmass
export acousticstiffness
export nzebcloadsacousticstiffness

include("FEMMAcoustSurfModule.jl")
using FinEtools.FEMMAcoustSurfModule
export FEMMAcoustSurf
export acousticABC
export pressure2resultantforce
export pressure2resultanttorque

include("FEMMDeforLinearBaseModule.jl")
using FinEtools.FEMMDeforLinearBaseModule
export FEMMDeforLinearAbstract
export stiffness
export nzebcloadsstiffness
export thermalstrainloads
export mass
export inspectintegpoints

include("FEMMDeforLinearModule.jl")
using FinEtools.FEMMDeforLinearModule
export FEMMDeforLinear
export stiffness
export nzebcloadsstiffness
export thermalstrainloads
export mass
export inspectintegpoints

include("FEMMDeforWinklerModule.jl")
using FinEtools.FEMMDeforWinklerModule
export FEMMDeforWinkler
export surfacenormalspringstiffness

include("FEMMDeforLinearMSModule.jl")
using FinEtools.FEMMDeforLinearMSModule
export FEMMDeforLinearMSH8, FEMMDeforLinearMST10
export stiffness
export nzebcloadsstiffness
export thermalstrainloads
export mass
export inspectintegpoints

include("AlgoBaseModule.jl")
using FinEtools.AlgoBaseModule
export FDataDict

include("AlgoAcoustModule.jl")
using FinEtools.AlgoAcoustModule

include("AlgoHeatDiffModule.jl")
using FinEtools.AlgoHeatDiffModule

include("AlgoDeforLinearModule.jl")
using FinEtools.AlgoDeforLinearModule

# println("")
# println("===================================================================")
# println("FinEtools ready")

end
