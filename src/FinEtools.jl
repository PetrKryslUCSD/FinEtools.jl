"""
FinEtools (C) 2017, Petr Krysl

Finite Element tools.  Julia implementation  of the finite element method
for continuum mechanics.
"""
module FinEtools

include("allmodules.jl")

###########################################################################
# General facilities
###########################################################################
using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
# Exported: basic numerical types, type of data dictionary
export FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict

using FinEtools.BoxModule: inbox, initbox!, updatebox!, boundingbox, inflatebox!, boxesoverlap, intersectboxes
# Exported: methods for manipulating and testing boxes
export inbox, initbox!, updatebox!, boundingbox, inflatebox!, boxesoverlap, intersectboxes

using FinEtools.PhysicalUnitModule: phun
# Exported: function for evaluating physical units
export phun

using FinEtools.RotationUtilModule: rotmat3!, skewmat!, cross3!, cross2
# Exported: functions for 3D rotation matrix computation, skew symmetric matrix computation,  in-place  cross product of 3-vectors, and  cross product of 2-vectors
export rotmat3!, skewmat!, cross3!, cross2

using FinEtools.CSysModule: CSys, updatecsmat!
# Exported: type  for coordinate systems, methods to invoke the update callback
export CSys, updatecsmat!

using FinEtools.FESetModule: FESet,  FESet0Manifold,  FESet1Manifold,  FESet2Manifold,  FESet3Manifold, manifdim, nodesperelem, count, fromarray!, connasarray, setlabel!, subset, cat, updateconn!, bfun, bfundpar, map2parametric, inparametric, centroidparametric,  FESetP1, FESetL2, FESetL3, FESetT3, FESetQ4, FESetQ9, FESetQ8, FESetT6, FESetH8, FESetH20, FESetH27, FESetT4, FESetT10
# Exported: types of finite element sets, abstract and concrete
export FESet,  FESet0Manifold,  FESet1Manifold,  FESet2Manifold,  FESet3Manifold, FESetP1, FESetL2, FESetL3, FESetT3, FESetQ4, FESetQ9, FESetQ8, FESetT6, FESetH8, FESetH20, FESetH27, FESetT4, FESetT10
# Exported: methods for accessing dimensions and counts
export manifdim, nodesperelem, count
# Exported: methods for  manipulating connectivity  and labels
export fromarray!, connasarray, setlabel!, subset, cat, updateconn!
# Exported: methods for computing basis function values and derivatives of basis functions  with respect to the parametric coordinates, and methods for working with parametric coordinates
export bfun, bfundpar, map2parametric, inparametric, centroidparametric 

using FinEtools.FENodeSetModule: FENodeSet, spacedim, xyz3, count
# Exported: type for FE node sets, methods for accessing dimensions  and counts
export FENodeSet, spacedim, xyz3, count

using FinEtools.FENodeToFEMapModule: FENodeToFEMap
# Exported: type for maps from nodes to finite elements
export FENodeToFEMap

using FinEtools.FieldModule: Field, ndofs,  nents, gathersysvec, gathersysvec!, gathervalues_asvec!, gathervalues_asmat!, gatherdofnums!, gatherfixedvalues_asvec!, gatherfixedvalues_asmat!,numberdofs!, setebc!, applyebc!, scattersysvec!, copyto!, wipe!
# Exported: abstract field type, methods for the abstract field type (retrieval of data  from a field, setting of data in the field)
export Field, ndofs,  nents, gathersysvec, gathersysvec!, gathervalues_asvec!, gathervalues_asmat!, gatherdofnums!, gatherfixedvalues_asvec!, gatherfixedvalues_asmat!, numberdofs!, setebc!, applyebc!, scattersysvec!, copyto!, wipe!

using FinEtools.GeneralFieldModule: GeneralField
# Exported: type of general field
export GeneralField

using FinEtools.NodalFieldModule: NodalField, nnodes
# Exported: type of nodal field
export NodalField, nnodes

using FinEtools.ElementalFieldModule: ElementalField, nelems
# Exported: type of elemental field
export ElementalField, nelems

using FinEtools.MeshUtilModule: linearspace, gradedspace
# Exported: functions to generate a sequence of numbers between start and stop
export linearspace, gradedspace

using FinEtools.MeshSelectionModule: connectednodes, connectedelems, selectnode, selectelem, findunconnnodes
# Exported: functions to select (find) nodes and elements
export connectednodes, connectedelems, selectnode, selectelem, findunconnnodes

using FinEtools.MeshExportModule: vtkexportmesh
# Exported: VTK export
export vtkexportmesh
using FinEtools.MeshExportModule: AbaqusExporter, close, HEADING, COMMENT, PART, END_PART,    ASSEMBLY, END_ASSEMBLY, INSTANCE, END_INSTANCE, NODE, ELEMENT,    NSET_NSET, ELSET_ELSET, ORIENTATION, MATERIAL, ELASTIC, EXPANSION, DENSITY, SECTION_CONTROLS, SOLID_SECTION, SURFACE_SECTION, STEP_PERTURBATION_STATIC, STEP_FREQUENCY,   STEP_PERTURBATION_BUCKLE, BOUNDARY, DLOAD, CLOAD, TEMPERATURE,    END_STEP,  NODE_PRINT, EL_PRINT,  ENERGY_PRINT
# Exported: Abaqus export
export AbaqusExporter, close, HEADING, COMMENT, PART, END_PART,    ASSEMBLY, END_ASSEMBLY, INSTANCE, END_INSTANCE, NODE, ELEMENT,    NSET_NSET, ELSET_ELSET, ORIENTATION, MATERIAL, ELASTIC, EXPANSION, DENSITY, SECTION_CONTROLS, SOLID_SECTION, SURFACE_SECTION, STEP_PERTURBATION_STATIC, STEP_FREQUENCY,   STEP_PERTURBATION_BUCKLE, BOUNDARY, DLOAD, CLOAD, TEMPERATURE,    END_STEP,  NODE_PRINT, EL_PRINT,  ENERGY_PRINT
using FinEtools.MeshExportModule: savecsv
# Exported: simple CSV export
export savecsv
using FinEtools.MeshExportModule: NASTRANExporter, close, CEND, BEGIN_BULK, ENDDATA, GRID, PSOLID, MAT1, CTETRA
# Exported: NASTRAN  export
export NASTRANExporter, close, CEND, BEGIN_BULK, ENDDATA, GRID, PSOLID, MAT1, CTETRA
using FinEtools.MeshExportModule: STLExporter, solid, facet, endsolid
# Exported: STL export
export STLExporter, solid, facet, endsolid
using FinEtools.MeshExportModule: h2libexporttri
# Exported: H2Lib export
export h2libexporttri

using FinEtools.MeshModificationModule: meshboundary,  fusenodes,  compactnodes,  mergemeshes, mergenmeshes, mergenodes,  renumberconn!,  meshsmoothing, mirrormesh, nodepartitioning, interior2boundary
# Exported: extraction of boundary, fusing of nodes  and merging of meshes, mesh smoothing,  node partitioning
export  meshboundary,  fusenodes,  compactnodes,  mergemeshes, mergenmeshes, mergenodes,  renumberconn!,  meshsmoothing, mirrormesh, nodepartitioning, interior2boundary

using FinEtools.MeshImportModule: import_NASTRAN, import_ABAQUS
# Exported: mesh import functions
export import_NASTRAN, import_ABAQUS

using FinEtools.VectorCacheModule: VectorCache, updateretrieve!
# Exported: vector-cache type and methods to invoke the update callback
export VectorCache, updateretrieve!

using FinEtools.ForceIntensityModule: ForceIntensity, updateforce!
# Exported: force-intensity type and methods to invoke the update callback
export ForceIntensity, updateforce!

using FinEtools.SurfaceNormalModule: SurfaceNormal, updatenormal!
# Exported: surface-normal evaluator type and methods to invoke the update callback
export SurfaceNormal, updatenormal!

using FinEtools.AssemblyModule: SysmatAssemblerBase, SysmatAssemblerSparse, SysmatAssemblerSparseSymm, SysmatAssemblerSparseHRZLumpingSymm, startassembly!, assemble!, makematrix!, SysvecAssemblerBase, SysvecAssembler, startassembly!, assemble!, makevector!
# Exported: types and methods for  sparse matrix assembly  and vector assembly
export SysmatAssemblerBase, SysmatAssemblerSparse, SysmatAssemblerSparseSymm, SysmatAssemblerSparseHRZLumpingSymm, startassembly!, assemble!, makematrix!, SysvecAssemblerBase, SysvecAssembler, startassembly!, assemble!, makevector!

using FinEtools.IntegRuleModule: IntegRule, TriRule, GaussRule, TetRule, PointRule, SimplexRule, TrapezoidalRule
# Exported: type for various integration rules
export IntegRule, TriRule, GaussRule, TetRule, PointRule, SimplexRule, TrapezoidalRule

using FinEtools.IntegDataModule: IntegData, otherdimensionunity, Jacobianpoint, Jacobiancurve, Jacobiansurface, Jacobianvolume, Jacobianmdim, integrationdata
# Exported: type to handle  integration data for various manifold dimensions
export IntegData, otherdimensionunity, Jacobianpoint, Jacobiancurve, Jacobiansurface, Jacobianvolume, Jacobianmdim, integrationdata

using FinEtools.FEMMBaseModule: FEMMAbstractBase, FEMMBase, associategeometry!, integratefieldfunction, integratefunction, transferfield!, distribloads, connectionmatrix, fieldfromintegpoints, elemfieldfromintegpoints, innerproduct
# Exported: type base discretization methods
export FEMMAbstractBase, FEMMBase, associategeometry!, integratefieldfunction, integratefunction, transferfield!, distribloads, connectionmatrix, fieldfromintegpoints, elemfieldfromintegpoints, innerproduct

###########################################################################
# Mesh-generation functionality for various shapes
###########################################################################
using FinEtools.MeshQuadrilateralModule: Q4annulus, Q4quadrilateral, Q4elliphole, Q4block, Q4blockx, Q4refine, Q8block, Q4toQ8, Q8annulus, Q8blockx
# Exported: mesh generation functions for quadrilateral elements
export Q4annulus, Q4quadrilateral, Q4elliphole, Q4block, Q4blockx, Q4refine, Q8block, Q4toQ8, Q8annulus, Q8blockx

using FinEtools.MeshLineModule: L2block, L2blockx, L3blockx
# Exported: mesh generation functions for line elements
export L2block, L2blockx, L3blockx

using FinEtools.MeshTriangleModule: T3blockx, T3block,  T3toT6,  T6block,  Q4toT3,  T3refine, T6blockx
# Exported: mesh generation functions for triangular elements
export  T3blockx, T3block,  T3toT6,  T6block,  Q4toT3,  T3refine, T6blockx

using FinEtools.MeshHexahedronModule: H8block,  H8blockx,  H8sphere,  H8refine, H8hexahedron, H8extrudeQ4, H8spheren, H8voximg,  H8layeredplatex, H8elliphole, H8toH27,  H27block, H20block,  H8toH20, H20blockx, H27blockx
# Exported: mesh generation functions for hexahedral elements
export  H8block,  H8blockx,  H8sphere,  H8refine, H8hexahedron, H8extrudeQ4, H8spheren, H8voximg,  H8layeredplatex, H8elliphole, H8toH27,  H27block, H20block,  H8toH20, H20blockx, H27blockx
    
using FinEtools.MeshTetrahedronModule: T4block, T4blockx, T4toT10, T10block, T10blockx, T10layeredplatex, T4meshedges, T4voximg
# Exported: mesh generation functions for tetrahedral elements
export  T4block, T4blockx, T4toT10, T10block, T10blockx, T10layeredplatex, T4meshedges, T4voximg

using FinEtools.VoxelBoxModule: VoxelBoxVolume, voxeldims, size, fillvolume!, fillsolid!,  intersectionop, unionop, complementop, differenceop,  solidsphere, solidhalfspace, solidbox, solidcylinder, trim, pad, threshold,  vtkexport
# Exported: type for voxel-box data structure, query methods
export VoxelBoxVolume, voxeldims, size
# Exported: methods to set voxel values to generate geometry
export fillvolume!, fillsolid!,  intersectionop, unionop, complementop,differenceop,  solidsphere, solidhalfspace, solidbox, solidcylinder
# Exported: methods for  manipulation and visualization  of voxel boxes
export trim, pad, threshold,  vtkexport
    
using FinEtools.VoxelTetMeshingModule: ElementSizeWeightFunction, ImageMesher, mesh!, volumes 
# Exported: type for the image mesher, type for control of element size gradation, method for generating  the mesh and queries
export ImageMesher, ElementSizeWeightFunction, mesh!, volumes 

###########################################################################
# Heat diffusion functionality
###########################################################################
using FinEtools.MatHeatDiffModule: MatHeatDiff
# Exported: type of heat-diffusion  material
export MatHeatDiff

using FinEtools.FEMMHeatDiffModule: FEMMHeatDiff, conductivity, nzebcloadsconductivity, energy
# Exported: type  for linear heat diffusion and discretization methods
export FEMMHeatDiff, conductivity, nzebcloadsconductivity, energy

using FinEtools.FEMMHeatDiffSurfModule: FEMMHeatDiffSurf, surfacetransfer, surfacetransferloads, nzebcsurfacetransferloads
# Exported: type  for linear heat diffusion boundary conditions and discretization methods
export FEMMHeatDiffSurf, surfacetransfer, surfacetransferloads, nzebcsurfacetransferloads

###########################################################################
# Acoustics functionality
###########################################################################
using FinEtools.MatAcoustFluidModule: MatAcoustFluid
# Exported: type of acoustic fluid material
export MatAcoustFluid

using FinEtools.FEMMAcoustModule: FEMMAcoust, acousticmass, nzebcloadsacousticmass,
acousticstiffness, nzebcloadsacousticstiffness
# Exported: type for linear acoustics  and discretization methods
export FEMMAcoust, acousticmass, nzebcloadsacousticmass,
    acousticstiffness, nzebcloadsacousticstiffness

using FinEtools.FEMMAcoustSurfModule: FEMMAcoustSurf, acousticABC, pressure2resultantforce, pressure2resultanttorque
# Exported: type for acoustic absorbing boundary condition  and  transformation matrices from pressure  to resultants
export FEMMAcoustSurf, acousticABC, pressure2resultantforce, pressure2resultanttorque

###########################################################################
# Linear deformation functionality
###########################################################################
using FinEtools.DeforModelRedModule: DeforModelRed, DeforModelRed1D, DeforModelRed2DStrain,    DeforModelRed2DStress, DeforModelRed2DAxisymm, DeforModelRed3D, nstressstrain, nthermstrain, stresscomponentmap, Blmat!
# Exported: types  for model reduction in stress analysis
export DeforModelRed, DeforModelRed1D, DeforModelRed2DStrain,    DeforModelRed2DStress, DeforModelRed2DAxisymm, DeforModelRed3D
# Exported: num stresses/strains,  number of thermal strains, and map of  the numbering of stress components
export nstressstrain, nthermstrain, stresscomponentmap 
# Exported: strain-displacement matrix  for all model-reduction types
export Blmat!

using FinEtools.MatDeforModule: MatDefor, strain2x2tto3v!, strain3vto2x2t!, strain3x3tto6v!, strain6vto3x3t!, strain9vto3x3t!, strain3x3tto9v!, strain9vto6v!, strain6vto9v!, stress2x2to3v!,  stress3vto2x2t!, stress3vto3x3t!, stress4vto3x3t!, stress6vto3x3t!, stress3x3tto6v!, stress9vto6v!,  stress6vto9v!, rotstressvec
# Exported: abstract type for  models of deformation,  conversion methods  for strain and stress, transformations  of strain and stress
export MatDefor, strain2x2tto3v!, strain3vto2x2t!, strain3x3tto6v!, strain6vto3x3t!, strain9vto3x3t!, strain3x3tto9v!, strain9vto6v!, strain6vto9v!, stress2x2to3v!,  stress3vto2x2t!, stress3vto3x3t!, stress4vto3x3t!, stress6vto3x3t!, stress3x3tto6v!, stress9vto6v!,  stress6vto9v!, rotstressvec

using FinEtools.MatDeforElastIsoModule: MatDeforElastIso
# Exported: type of  isotropic elastic material
export MatDeforElastIso

using FinEtools.MatDeforElastOrthoModule: MatDeforElastOrtho
# Exported: type of orthotropic elastic material
export MatDeforElastOrtho

using FinEtools.FEMMDeforLinearBaseModule: FEMMDeforLinearAbstract, stiffness, nzebcloadsstiffness, thermalstrainloads, mass, inspectintegpoints
# Exported: abstract type for linear information, discretization methods for the abstract type 
export FEMMDeforLinearAbstract, stiffness, nzebcloadsstiffness, thermalstrainloads, mass, inspectintegpoints

using FinEtools.FEMMDeforLinearModule: FEMMDeforLinear
# Exported: type for linear deformation
export FEMMDeforLinear

using FinEtools.FEMMDeforWinklerModule: FEMMDeforWinkler, surfacenormalspringstiffness
# Exported: type for distributed-spring support, discretization method 
export FEMMDeforWinkler, surfacenormalspringstiffness

using FinEtools.FEMMDeforLinearMSModule: FEMMDeforLinearMSH8, FEMMDeforLinearMST10, stiffness, nzebcloadsstiffness, thermalstrainloads, inspectintegpoints
# Exported: type for mean-strain solid elements, discretization methods
export FEMMDeforLinearMSH8, FEMMDeforLinearMST10, stiffness, nzebcloadsstiffness, thermalstrainloads, inspectintegpoints

using FinEtools.FEMMDeforSurfaceDampingModule: FEMMDeforSurfaceDamping, dampingABC
#Exported: type for surface damping (absorbing boundary conditions)
export FEMMDeforSurfaceDamping, dampingABC

end
