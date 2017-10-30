"""
FinEtools (C) 2017, Petr Krysl

Finite Element tools.  Julia implementation  of the finite element method
for continuum mechanics.
"""
module FinEtools

include("FTypesModule.jl")
using FinEtools.FTypesModule
export FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec
export FDataDict

include("BoxModule.jl")
using FinEtools.BoxModule
export inbox, initbox!, updatebox!, boundingbox, inflatebox!, boxesoverlap

include("MatrixUtilityModule.jl")

include("PhysicalUnitModule.jl")
using FinEtools.PhysicalUnitModule
phun = FinEtools.PhysicalUnitModule.phun
export phun

include("RotationUtilModule.jl")
using FinEtools.RotationUtilModule
export rotmat,  rotmat3!, skewmat!, cross3, cross3!, cross2

include("CSysModule.jl")
using FinEtools.CSysModule
export CSys, updatecsmat!

include("FESetModule.jl")
using FinEtools.FESetModule
export FESet,  FESet0Manifold,  FESet1Manifold,  FESet2Manifold,  FESet3Manifold
export manifdim, nodesperelem, count, getconn!, setlabel!, subset, cat, updateconn!
export bfun, bfundpar, map2parametric, inparametric, centroidparametric
export FESetP1
export FESetL2, FESetL3
export FESetT3, FESetQ4, FESetQ9, FESetQ8, FESetT6
export FESetH8, FESetH20, FESetH27, FESetT4, FESetT10

include("FENodeSetModule.jl")
using FinEtools.FENodeSetModule
export FENodeSet, spacedim, xyz3, count

include("FENodeToFEMapModule.jl")
using FinEtools.FENodeToFEMapModule
export FENodeToFEMap

include("FieldModule.jl")
using FinEtools.FieldModule
export Field, ndofs,  nents, gathersysvec, gathersysvec!, gathervalues_asvec!,
    gathervalues_asmat!, gatherdofnums!, gatherfixedvalues_asvec!,
    gatherfixedvalues_asmat!, numberdofs!, setebc!, applyebc!,
    scattersysvec!, copy!, wipe!

include("GeneralFieldModule.jl")
using FinEtools.GeneralFieldModule
export GeneralField

include("NodalFieldModule.jl")
using FinEtools.NodalFieldModule
export NodalField, nnodes

include("ElementalFieldModule.jl")
using FinEtools.ElementalFieldModule
export ElementalField, nelems

include("MeshUtilModule.jl")
using FinEtools.MeshUtilModule

include("MeshSelectionModule.jl")
using FinEtools.MeshSelectionModule
export connectedelems, connectednodes, selectnode, selectelem, findunconnnodes

include("MeshExportModule.jl")
using FinEtools.MeshExportModule
export vtkexportmesh
export AbaqusExporter, close, HEADING, COMMENT, PART, END_PART,
    ASSEMBLY, END_ASSEMBLY, INSTANCE, END_INSTANCE, NODE, ELEMENT,
    NSET_NSET, ELSET_ELSET, ORIENTATION, MATERIAL,
    ELASTIC, EXPANSION, DENSITY, SECTION_CONTROLS, SOLID_SECTION,
    SURFACE_SECTION, STEP_PERTURBATION_STATIC, STEP_FREQUENCY,
    STEP_PERTURBATION_BUCKLE, BOUNDARY, DLOAD, CLOAD, TEMPERATURE,
    END_STEP,  NODE_PRINT, EL_PRINT,  ENERGY_PRINT
export savecsv
export export_NASTRAN

include("MeshImportModule.jl")
using FinEtools.MeshImportModule

include("MeshModificationModule.jl")
using FinEtools.MeshModificationModule
export  meshboundary,  fusenodes,  compactnodes,  mergemeshes,  mergenmeshes,
    mergenodes,  renumberconn!,  meshsmoothing,  mirrormesh, nodepartitioning, 
    interior2boundary

include("MeshQuadrilateralModule.jl")
using FinEtools.MeshQuadrilateralModule
export Q4annulus, Q4quadrilateral, Q4elliphole, Q4block, Q4blockx, Q4refine
export Q8block, Q4toQ8, Q8annulus, Q8blockx

include("MeshLineModule.jl")
using FinEtools.MeshLineModule
export L2block, L2blockx, L3blockx

include("MeshTriangleModule.jl")
using FinEtools.MeshTriangleModule
export  T3blockx, T3block,  T3toT6,  T6block,  Q4toT3,  T3refine, T6blockx

include("MeshHexahedronModule.jl")
using FinEtools.MeshHexahedronModule
export  H8block,  H8blockx,  H8sphere,  H8refine, H8hexahedron, H8extrudeQ4,
    H8spheren, H8voximg,  H8compositeplatex, H8elliphole, H8toH27,  H27block,
    H20block,  H8toH20, H20blockx, H27blockx
    
include("MeshTetrahedronModule.jl")
using FinEtools.MeshTetrahedronModule
export  T4block, T4blockx, T4toT10, T10block, T10blockx, T10compositeplatex, T4meshedges, T4voximg

include("ForceIntensityModule.jl")
using FinEtools.ForceIntensityModule
export ForceIntensity, getforce!

include("MatHeatDiffModule.jl")
using FinEtools.MatHeatDiffModule
export MatHeatDiff

include("MatAcoustFluidModule.jl")
using FinEtools.MatAcoustFluidModule
export MatAcoustFluid

include("DeforModelRedModule.jl")
using FinEtools.DeforModelRedModule
export DeforModelRed, DeforModelRed1D, DeforModelRed2DStrain,
    DeforModelRed2DStress, DeforModelRed2DAxisymm, DeforModelRed3D
export nstrains, stresscomponentmap, Blmat!

include("MatDeforModule.jl")
using FinEtools.MatDeforModule
export MatDefor
export strain2x2tto3v!, strain3vto2x2t!, strain3x3tto6v!, strain6vto3x3t!,
    strain9vto3x3t!, strain3x3tto9v!, strain9vto6v!, strain6vto9v!
export stress2x2to3v!,  stress3vto2x2t!, stress3vto3x3t!, stress4vto3x3t!,
    stress6vto3x3t!, stress3x3tto6v!, stress9vto6v!,  stress6vto9v!
export rotstressvec

include("MatDeforElastIsoModule.jl")
using FinEtools.MatDeforElastIsoModule
export MatDeforElastIso

include("MatDeforElastOrthoModule.jl")
using FinEtools.MatDeforElastOrthoModule
export MatDeforElastOrtho

include("AssemblyModule.jl")
using FinEtools.AssemblyModule
export SysmatAssemblerBase, SysmatAssemblerSparse, SysmatAssemblerSparseSymm
export startassembly!, assemble!, makematrix!
export SysvecAssemblerBase, SysvecAssembler
export startassembly!, assemble!, makevector!

include("IntegRuleModule.jl")
using FinEtools.IntegRuleModule
export IntegRule, TriRule, GaussRule, TetRule, PointRule, SimplexRule

include("GeoDModule.jl")
using FinEtools.GeoDModule
export GeoD
export otherdimensionunity
export Jacobianpoint, Jacobiancurve, Jacobiansurface, Jacobianvolume, Jacobianmdim

include("FEMMBaseModule.jl")
using FinEtools.FEMMBaseModule
export FEMMAbstractBase, FEMMBase
export associategeometry!, integratefieldfunction, integratefunction,
    transferfield!, distribloads, connectionmatrix,
    fieldfromintegpoints, elemfieldfromintegpoints

include("FEMMHeatDiffModule.jl")
using FinEtools.FEMMHeatDiffModule
export FEMMHeatDiff
export conductivity, nzebcloadsconductivity

include("FEMMHeatDiffSurfModule.jl")
using FinEtools.FEMMHeatDiffSurfModule
export FEMMHeatDiffSurf
export surfacetransfer, surfacetransferloads, nzebcsurfacetransferloads

include("FEMMAcoustModule.jl")
using FinEtools.FEMMAcoustModule
export FEMMAcoust
export acousticmass, nzebcloadsacousticmass,
    acousticstiffness, nzebcloadsacousticstiffness

include("FEMMAcoustSurfModule.jl")
using FinEtools.FEMMAcoustSurfModule
export FEMMAcoustSurf
export acousticABC, pressure2resultantforce, pressure2resultanttorque

include("FEMMDeforLinearBaseModule.jl")
using FinEtools.FEMMDeforLinearBaseModule
export FEMMDeforLinearAbstract
export stiffness, nzebcloadsstiffness, thermalstrainloads, mass,
    inspectintegpoints

include("FEMMDeforLinearModule.jl")
using FinEtools.FEMMDeforLinearModule
export FEMMDeforLinear
export stiffness, nzebcloadsstiffness, thermalstrainloads, inspectintegpoints

include("FEMMDeforWinklerModule.jl")
using FinEtools.FEMMDeforWinklerModule
export FEMMDeforWinkler
export surfacenormalspringstiffness

include("FEMMDeforLinearMSModule.jl")
using FinEtools.FEMMDeforLinearMSModule
export FEMMDeforLinearMSH8, FEMMDeforLinearMST10
export stiffness, nzebcloadsstiffness, thermalstrainloads, inspectintegpoints

include("AlgoBaseModule.jl")
using FinEtools.AlgoBaseModule

include("AlgoAcoustModule.jl")
using FinEtools.AlgoAcoustModule

include("AlgoHeatDiffModule.jl")
using FinEtools.AlgoHeatDiffModule

include("AlgoDeforLinearModule.jl")
using FinEtools.AlgoDeforLinearModule


include("VoxelBoxModule.jl")
using FinEtools.VoxelBoxModule
export VoxelBoxVolume, voxeldims, size,
    fillvolume!, fillsolid!,
    intersectionop, unionop, complementop, differenceop,
    solidsphere, solidhalfspace, solidbox, solidcylinder,
    trim, pad, threshold,
    vtkexport
    
include("TetRemeshingModule.jl")
using FinEtools.TetRemeshingModule

include("VoxelTetMeshingModule.jl")
export ElementSizeWeightFunction, ImageMesher, mesh!, volumes 

end
