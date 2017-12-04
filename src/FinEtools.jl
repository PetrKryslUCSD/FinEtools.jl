"""
FinEtools (C) 2017, Petr Krysl

Finite Element tools.  Julia implementation  of the finite element method
for continuum mechanics.
"""
module FinEtools

include("FTypesModule.jl")
using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
export FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict

include("BoxModule.jl")
using FinEtools.BoxModule: inbox, initbox!, updatebox!, boundingbox, inflatebox!, boxesoverlap
export inbox, initbox!, updatebox!, boundingbox, inflatebox!, boxesoverlap

include("MatrixUtilityModule.jl")

include("PhysicalUnitModule.jl")
using FinEtools.PhysicalUnitModule: phun
export phun

include("RotationUtilModule.jl")
using FinEtools.RotationUtilModule: rotmat3!, skewmat!, cross3, cross3!, cross2
export rotmat3!, skewmat!, cross3, cross3!, cross2

include("CSysModule.jl")
using FinEtools.CSysModule: CSys, updatecsmat!
export CSys, updatecsmat!

include("FESetModule.jl")
using FinEtools.FESetModule: FESet,  FESet0Manifold,  FESet1Manifold,  FESet2Manifold,  FESet3Manifold, manifdim, nodesperelem, count, fromarray!, connasarray, setlabel!, subset, cat, updateconn!, bfun, bfundpar, map2parametric, inparametric, centroidparametric,  FESetP1, FESetL2, FESetL3, FESetT3, FESetQ4, FESetQ9, FESetQ8, FESetT6, FESetH8, FESetH20, FESetH27, FESetT4, FESetT10
export FESet,  FESet0Manifold,  FESet1Manifold,  FESet2Manifold,  FESet3Manifold, manifdim, nodesperelem, count, fromarray!, connasarray, setlabel!, subset, cat, updateconn!, bfun, bfundpar, map2parametric, inparametric, centroidparametric,  FESetP1, FESetL2, FESetL3, FESetT3, FESetQ4, FESetQ9, FESetQ8, FESetT6, FESetH8, FESetH20, FESetH27, FESetT4, FESetT10

include("FENodeSetModule.jl")
using FinEtools.FENodeSetModule: FENodeSet, spacedim, xyz3, count
export FENodeSet, spacedim, xyz3, count

include("FENodeToFEMapModule.jl")
using FinEtools.FENodeToFEMapModule: FENodeToFEMap
export FENodeToFEMap

include("FieldModule.jl")
using FinEtools.FieldModule: Field, ndofs,  nents, gathersysvec, gathersysvec!, gathervalues_asvec!, gathervalues_asmat!, gatherdofnums!, gatherfixedvalues_asvec!, gatherfixedvalues_asmat!,numberdofs!, setebc!, applyebc!, scattersysvec!, copy!, wipe!
export Field, ndofs,  nents, gathersysvec, gathersysvec!, gathervalues_asvec!, gathervalues_asmat!, gatherdofnums!, gatherfixedvalues_asvec!, gatherfixedvalues_asmat!, numberdofs!, setebc!, applyebc!, scattersysvec!, copy!, wipe!

include("GeneralFieldModule.jl")
using FinEtools.GeneralFieldModule: GeneralField
export GeneralField

include("NodalFieldModule.jl")
using FinEtools.NodalFieldModule: NodalField, nnodes
export NodalField, nnodes

include("ElementalFieldModule.jl")
using FinEtools.ElementalFieldModule: ElementalField, nelems
export ElementalField, nelems

include("MeshUtilModule.jl")
using FinEtools.MeshUtilModule

include("MeshSelectionModule.jl")
using FinEtools.MeshSelectionModule: connectednodes, selectnode, selectelem, findunconnnodes
export connectednodes, selectnode, selectelem, findunconnnodes

include("MeshExportModule.jl")
using FinEtools.MeshExportModule: vtkexportmesh
export vtkexportmesh
using FinEtools.MeshExportModule: AbaqusExporter, close, HEADING, COMMENT, PART, END_PART,    ASSEMBLY, END_ASSEMBLY, INSTANCE, END_INSTANCE, NODE, ELEMENT,    NSET_NSET, ELSET_ELSET, ORIENTATION, MATERIAL, ELASTIC, EXPANSION, DENSITY, SECTION_CONTROLS, SOLID_SECTION, SURFACE_SECTION, STEP_PERTURBATION_STATIC, STEP_FREQUENCY,   STEP_PERTURBATION_BUCKLE, BOUNDARY, DLOAD, CLOAD, TEMPERATURE,    END_STEP,  NODE_PRINT, EL_PRINT,  ENERGY_PRINT
export AbaqusExporter, close, HEADING, COMMENT, PART, END_PART,    ASSEMBLY, END_ASSEMBLY, INSTANCE, END_INSTANCE, NODE, ELEMENT,    NSET_NSET, ELSET_ELSET, ORIENTATION, MATERIAL, ELASTIC, EXPANSION, DENSITY, SECTION_CONTROLS, SOLID_SECTION, SURFACE_SECTION, STEP_PERTURBATION_STATIC, STEP_FREQUENCY,   STEP_PERTURBATION_BUCKLE, BOUNDARY, DLOAD, CLOAD, TEMPERATURE,    END_STEP,  NODE_PRINT, EL_PRINT,  ENERGY_PRINT
using FinEtools.MeshExportModule: savecsv
export savecsv
using FinEtools.MeshExportModule: NASTRANExporter, close, CEND, BEGIN_BULK, ENDDATA, GRID, PSOLID, MAT1, CTETRA
export NASTRANExporter, close, CEND, BEGIN_BULK, ENDDATA, GRID, PSOLID, MAT1, CTETRA
using FinEtools.MeshExportModule: STLExporter, solid, facet, endsolid
export STLExporter, solid, facet, endsolid

include("MeshImportModule.jl")
using FinEtools.MeshImportModule

include("MeshModificationModule.jl")
using FinEtools.MeshModificationModule: meshboundary,  fusenodes,  compactnodes,  mergemeshes, mergenmeshes, mergenodes,  renumberconn!,  meshsmoothing, mirrormesh, nodepartitioning, interior2boundary
export  meshboundary,  fusenodes,  compactnodes,  mergemeshes, mergenmeshes, mergenodes,  renumberconn!,  meshsmoothing, mirrormesh, nodepartitioning, interior2boundary

include("MeshQuadrilateralModule.jl")
using FinEtools.MeshQuadrilateralModule: Q4annulus, Q4quadrilateral, Q4elliphole, Q4block, Q4blockx, Q4refine, Q8block, Q4toQ8, Q8annulus, Q8blockx
export Q4annulus, Q4quadrilateral, Q4elliphole, Q4block, Q4blockx, Q4refine, Q8block, Q4toQ8, Q8annulus, Q8blockx

include("MeshLineModule.jl")
using FinEtools.MeshLineModule: L2block, L2blockx, L3blockx
export L2block, L2blockx, L3blockx

include("MeshTriangleModule.jl")
using FinEtools.MeshTriangleModule: T3blockx, T3block,  T3toT6,  T6block,  Q4toT3,  T3refine, T6blockx
export  T3blockx, T3block,  T3toT6,  T6block,  Q4toT3,  T3refine, T6blockx

include("MeshHexahedronModule.jl")
using FinEtools.MeshHexahedronModule: H8block,  H8blockx,  H8sphere,  H8refine, H8hexahedron, H8extrudeQ4, H8spheren, H8voximg,  H8layeredplatex, H8elliphole, H8toH27,  H27block, H20block,  H8toH20, H20blockx, H27blockx
export  H8block,  H8blockx,  H8sphere,  H8refine, H8hexahedron, H8extrudeQ4, H8spheren, H8voximg,  H8layeredplatex, H8elliphole, H8toH27,  H27block, H20block,  H8toH20, H20blockx, H27blockx
    
include("MeshTetrahedronModule.jl")
using FinEtools.MeshTetrahedronModule: T4block, T4blockx, T4toT10, T10block, T10blockx, T10layeredplatex, T4meshedges, T4voximg
export  T4block, T4blockx, T4toT10, T10block, T10blockx, T10layeredplatex, T4meshedges, T4voximg

include("ForceIntensityModule.jl")
using FinEtools.ForceIntensityModule: ForceIntensity, updateforce!
export ForceIntensity, updateforce!

include("MatHeatDiffModule.jl")
using FinEtools.MatHeatDiffModule: MatHeatDiff
export MatHeatDiff

include("MatAcoustFluidModule.jl")
using FinEtools.MatAcoustFluidModule: MatAcoustFluid
export MatAcoustFluid

include("DeforModelRedModule.jl")
using FinEtools.DeforModelRedModule: DeforModelRed, DeforModelRed1D, DeforModelRed2DStrain,    DeforModelRed2DStress, DeforModelRed2DAxisymm, DeforModelRed3D, nstressstrain, nthermstrain, stresscomponentmap, Blmat!
export DeforModelRed, DeforModelRed1D, DeforModelRed2DStrain,    DeforModelRed2DStress, DeforModelRed2DAxisymm, DeforModelRed3D, nstressstrain, nthermstrain, stresscomponentmap, Blmat!

include("MatDeforModule.jl")
using FinEtools.MatDeforModule: MatDefor, strain2x2tto3v!, strain3vto2x2t!, strain3x3tto6v!, strain6vto3x3t!, strain9vto3x3t!, strain3x3tto9v!, strain9vto6v!, strain6vto9v!, stress2x2to3v!,  stress3vto2x2t!, stress3vto3x3t!, stress4vto3x3t!, stress6vto3x3t!, stress3x3tto6v!, stress9vto6v!,  stress6vto9v!, rotstressvec
export MatDefor, strain2x2tto3v!, strain3vto2x2t!, strain3x3tto6v!, strain6vto3x3t!, strain9vto3x3t!, strain3x3tto9v!, strain9vto6v!, strain6vto9v!, stress2x2to3v!,  stress3vto2x2t!, stress3vto3x3t!, stress4vto3x3t!, stress6vto3x3t!, stress3x3tto6v!, stress9vto6v!,  stress6vto9v!, rotstressvec

include("MatDeforElastIsoModule.jl")
using FinEtools.MatDeforElastIsoModule: MatDeforElastIso
export MatDeforElastIso

include("MatDeforElastOrthoModule.jl")
using FinEtools.MatDeforElastOrthoModule: MatDeforElastOrtho
export MatDeforElastOrtho

include("AssemblyModule.jl")
using FinEtools.AssemblyModule: SysmatAssemblerBase, SysmatAssemblerSparse, SysmatAssemblerSparseSymm, startassembly!, assemble!, makematrix!, SysvecAssemblerBase, SysvecAssembler, startassembly!, assemble!, makevector!
export SysmatAssemblerBase, SysmatAssemblerSparse, SysmatAssemblerSparseSymm, startassembly!, assemble!, makematrix!, SysvecAssemblerBase, SysvecAssembler, startassembly!, assemble!, makevector!

include("IntegRuleModule.jl")
using FinEtools.IntegRuleModule: IntegRule, TriRule, GaussRule, TetRule, PointRule, SimplexRule
export IntegRule, TriRule, GaussRule, TetRule, PointRule, SimplexRule

include("IntegDataModule.jl")
using FinEtools.IntegDataModule: IntegData, otherdimensionunity, Jacobianpoint, Jacobiancurve, Jacobiansurface, Jacobianvolume, Jacobianmdim, integrationdata
export IntegData, otherdimensionunity, Jacobianpoint, Jacobiancurve, Jacobiansurface, Jacobianvolume, Jacobianmdim, integrationdata

include("FEMMBaseModule.jl")
using FinEtools.FEMMBaseModule: FEMMAbstractBase, FEMMBase, associategeometry!, integratefieldfunction, integratefunction, transferfield!, distribloads, connectionmatrix, fieldfromintegpoints, elemfieldfromintegpoints
export FEMMAbstractBase, FEMMBase, associategeometry!, integratefieldfunction, integratefunction, transferfield!, distribloads, connectionmatrix, fieldfromintegpoints, elemfieldfromintegpoints

include("FEMMHeatDiffModule.jl")
using FinEtools.FEMMHeatDiffModule: FEMMHeatDiff, conductivity, nzebcloadsconductivity
export FEMMHeatDiff, conductivity, nzebcloadsconductivity

include("FEMMHeatDiffSurfModule.jl")
using FinEtools.FEMMHeatDiffSurfModule: FEMMHeatDiffSurf, surfacetransfer, surfacetransferloads, nzebcsurfacetransferloads
export FEMMHeatDiffSurf, surfacetransfer, surfacetransferloads, nzebcsurfacetransferloads

include("FEMMAcoustModule.jl")
using FinEtools.FEMMAcoustModule: FEMMAcoust, acousticmass, nzebcloadsacousticmass,
acousticstiffness, nzebcloadsacousticstiffness
export FEMMAcoust, acousticmass, nzebcloadsacousticmass,
    acousticstiffness, nzebcloadsacousticstiffness

include("FEMMAcoustSurfModule.jl")
using FinEtools.FEMMAcoustSurfModule: FEMMAcoustSurf, acousticABC, pressure2resultantforce, pressure2resultanttorque
export FEMMAcoustSurf, acousticABC, pressure2resultantforce, pressure2resultanttorque

include("FEMMDeforLinearBaseModule.jl")
using FinEtools.FEMMDeforLinearBaseModule: FEMMDeforLinearAbstract, stiffness, nzebcloadsstiffness, thermalstrainloads, mass, inspectintegpoints
export FEMMDeforLinearAbstract, stiffness, nzebcloadsstiffness, thermalstrainloads, mass, inspectintegpoints

include("FEMMDeforLinearModule.jl")
using FinEtools.FEMMDeforLinearModule: FEMMDeforLinear, stiffness, nzebcloadsstiffness, thermalstrainloads, inspectintegpoints
export FEMMDeforLinear, stiffness, nzebcloadsstiffness, thermalstrainloads, inspectintegpoints

include("FEMMDeforWinklerModule.jl")
using FinEtools.FEMMDeforWinklerModule: FEMMDeforWinkler, surfacenormalspringstiffness
export FEMMDeforWinkler, surfacenormalspringstiffness

include("FEMMDeforLinearMSModule.jl")
using FinEtools.FEMMDeforLinearMSModule: FEMMDeforLinearMSH8, FEMMDeforLinearMST10, stiffness, nzebcloadsstiffness, thermalstrainloads, inspectintegpoints
export FEMMDeforLinearMSH8, FEMMDeforLinearMST10, stiffness, nzebcloadsstiffness, thermalstrainloads, inspectintegpoints

include("AlgoBaseModule.jl")
using FinEtools.AlgoBaseModule

include("AlgoAcoustModule.jl")
using FinEtools.AlgoAcoustModule

include("AlgoHeatDiffModule.jl")
using FinEtools.AlgoHeatDiffModule

include("AlgoDeforLinearModule.jl")
using FinEtools.AlgoDeforLinearModule


include("VoxelBoxModule.jl")
using FinEtools.VoxelBoxModule: VoxelBoxVolume, voxeldims, size, fillvolume!, fillsolid!,  intersectionop, unionop, complementop, differenceop,  solidsphere, solidhalfspace, solidbox, solidcylinder, trim, pad, threshold,  vtkexport
export VoxelBoxVolume, voxeldims, size, fillvolume!, fillsolid!,  intersectionop, unionop, complementop, differenceop,  solidsphere, solidhalfspace, solidbox, solidcylinder, trim, pad, threshold,  vtkexport
    
include("TetRemeshingModule.jl")
using FinEtools.TetRemeshingModule

include("VoxelTetMeshingModule.jl")
using FinEtools.VoxelTetMeshingModule: ElementSizeWeightFunction, ImageMesher, mesh!, volumes 
export ElementSizeWeightFunction, ImageMesher, mesh!, volumes 

end
