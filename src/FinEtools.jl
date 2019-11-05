"""
FinEtools (C) 2017-2019, Petr Krysl

Finite Element tools.  Julia implementation  of the finite element method
for continuum mechanics.
"""
module FinEtools

__precompile__(true)

include("allmodules.jl")

# Exports follow:

###########################################################################
# General facilities
###########################################################################
using .FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
# Exported: basic numerical types, type of data dictionary
export FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict

using .BoxModule: inbox, initbox!, updatebox!, boundingbox, inflatebox!, boxesoverlap, intersectboxes
# Exported: methods for manipulating and testing boxes
export inbox, initbox!, updatebox!, boundingbox, inflatebox!, boxesoverlap, intersectboxes

using .PhysicalUnitModule: phun
# Exported: function for evaluating physical units
export phun

using .RotationUtilModule: rotmat3!, skewmat!, cross3!, cross2
# Exported: functions for 3D rotation matrix computation, skew symmetric matrix computation,  in-place  cross product of 3-vectors, and  cross product of 2-vectors
export rotmat3!, skewmat!, cross3!, cross2

using .CSysModule: CSys, updatecsmat!
# Exported: type  for coordinate systems, methods to invoke the update callback
export CSys, updatecsmat!

using .FESetModule: AbstractFESet,  AbstractFESet0Manifold,  AbstractFESet1Manifold,  AbstractFESet2Manifold,  AbstractFESet3Manifold, manifdim, nodesperelem, count, fromarray!, connasarray, setlabel!, subset, cat, updateconn!, bfun, bfundpar, map2parametric, inparametric, centroidparametric,  FESetP1, FESetL2, FESetL3, FESetT3, FESetQ4, FESetQ9, FESetQ8, FESetT6, FESetH8, FESetH20, FESetH27, FESetT4, FESetT10
# Exported: types of finite element sets, abstract and concrete
export AbstractFESet,  AbstractFESet0Manifold,  AbstractFESet1Manifold,  AbstractFESet2Manifold,  AbstractFESet3Manifold, FESetP1, FESetL2, FESetL3, FESetT3, FESetQ4, FESetQ9, FESetQ8, FESetT6, FESetH8, FESetH20, FESetH27, FESetT4, FESetT10
# Exported: methods for accessing dimensions and counts
export manifdim, nodesperelem, count
# Exported: methods for  manipulating connectivity  and labels
export fromarray!, connasarray, setlabel!, subset, cat, updateconn!
# Exported: methods for computing basis function values and derivatives of basis functions  with respect to the parametric coordinates, and methods for working with parametric coordinates
export bfun, bfundpar, map2parametric, inparametric, centroidparametric

using .FENodeSetModule: FENodeSet, spacedim, xyz3, count
# Exported: type for FE node sets, methods for accessing dimensions  and counts
export FENodeSet, spacedim, xyz3, count

using .FENodeToFEMapModule: FENodeToFEMap
# Exported: type for maps from nodes to finite elements
export FENodeToFEMap

using .FieldModule: AbstractField, ndofs,  nents, gathersysvec, gathersysvec!, gathervalues_asvec!, gathervalues_asmat!, gatherdofnums!, gatherfixedvalues_asvec!, gatherfixedvalues_asmat!,numberdofs!, setebc!, applyebc!, scattersysvec!, copyto!, wipe!, prescribeddofs
# Exported: abstract field type, methods for the abstract field type (retrieval of data  from a field, setting of data in the field)
export AbstractField, ndofs,  nents, gathersysvec, gathersysvec!, gathervalues_asvec!, gathervalues_asmat!, gatherdofnums!, gatherfixedvalues_asvec!, gatherfixedvalues_asmat!, numberdofs!, setebc!, applyebc!, scattersysvec!, copyto!, wipe!, prescribeddofs

using .GeneralFieldModule: GeneralField
# Exported: type of general field
export GeneralField

using .NodalFieldModule: NodalField, nnodes
# Exported: type of nodal field
export NodalField, nnodes

using .ElementalFieldModule: ElementalField, nelems
# Exported: type of elemental field
export ElementalField, nelems

using .MeshUtilModule: linearspace, gradedspace
# Exported: functions to generate a sequence of numbers between start and stop
export linearspace, gradedspace

using .MeshSelectionModule: connectednodes, connectedelems, selectnode, selectelem, findunconnnodes
# Exported: functions to select (find) nodes and elements
export connectednodes, connectedelems, selectnode, selectelem, findunconnnodes

using .MeshExportModule.VTK: vtkexportmesh
# Exported: VTK export
export vtkexportmesh
using .MeshExportModule.Abaqus: AbaqusExporter, close, HEADING, COMMENT, PART, END_PART,    ASSEMBLY, END_ASSEMBLY, INSTANCE, END_INSTANCE, NODE, ELEMENT,    NSET_NSET, ELSET_ELSET, ORIENTATION, MATERIAL, ELASTIC, EXPANSION, DENSITY, SECTION_CONTROLS, SOLID_SECTION, SURFACE_SECTION, STEP_PERTURBATION_STATIC, STEP_FREQUENCY,   STEP_PERTURBATION_BUCKLE, BOUNDARY, DLOAD, CLOAD, TEMPERATURE,    END_STEP,  NODE_PRINT, EL_PRINT,  ENERGY_PRINT
# Exported: Abaqus export
export AbaqusExporter, close, HEADING, COMMENT, PART, END_PART,    ASSEMBLY, END_ASSEMBLY, INSTANCE, END_INSTANCE, NODE, ELEMENT,    NSET_NSET, ELSET_ELSET, ORIENTATION, MATERIAL, ELASTIC, EXPANSION, DENSITY, SECTION_CONTROLS, SOLID_SECTION, SURFACE_SECTION, STEP_PERTURBATION_STATIC, STEP_FREQUENCY,   STEP_PERTURBATION_BUCKLE, BOUNDARY, DLOAD, CLOAD, TEMPERATURE,    END_STEP,  NODE_PRINT, EL_PRINT,  ENERGY_PRINT
using .MeshExportModule.CSV: savecsv
# Exported: simple CSV export
export savecsv
using .MeshExportModule.NASTRAN: NASTRANExporter, close, CEND, BEGIN_BULK, ENDDATA, GRID, PSOLID, MAT1, CTETRA
# Exported: NASTRAN  export
export NASTRANExporter, close, CEND, BEGIN_BULK, ENDDATA, GRID, PSOLID, MAT1, CTETRA
using .MeshExportModule.STL: STLExporter, solid, facet, endsolid
# Exported: STL export
export STLExporter, solid, facet, endsolid
using .MeshExportModule.H2Lib: h2libexporttri
# Exported: H2Lib export
export h2libexporttri

using .MeshModificationModule: meshboundary,  fusenodes,  compactnodes,  mergemeshes, mergenmeshes, mergenodes,  renumberconn!,  meshsmoothing, mirrormesh, nodepartitioning, interior2boundary, adjgraph, nodedegrees, revcm
# Exported: extraction of boundary, fusing of nodes  and merging of meshes, mesh smoothing,  node partitioning
export  meshboundary,  fusenodes,  compactnodes,  mergemeshes, mergenmeshes, mergenodes,  renumberconn!,  meshsmoothing, mirrormesh, nodepartitioning, interior2boundary, adjgraph, nodedegrees, revcm

using .MeshImportModule: import_NASTRAN, import_ABAQUS
# Exported: mesh import functions
export import_NASTRAN, import_ABAQUS

using .VectorCacheModule: VectorCache, updateretrieve!, settime!
# Exported: vector-cache type and methods to invoke the update callback
export VectorCache, updateretrieve!, settime!

using .ForceIntensityModule: ForceIntensity, updateforce!, settime!
# Exported: force-intensity type and methods to invoke the update callback
export ForceIntensity, updateforce!, settime!

using .SurfaceNormalModule: SurfaceNormal, updatenormal!
# Exported: surface-normal evaluator type and methods to invoke the update callback
export SurfaceNormal, updatenormal!

using .AssemblyModule: AbstractSysmatAssembler, SysmatAssemblerSparse, SysmatAssemblerSparseSymm, SysmatAssemblerSparseHRZLumpingSymm, startassembly!, assemble!, makematrix!, AbstractSysvecAssembler, SysvecAssembler, startassembly!, assemble!, makevector!
# Exported: types and methods for  sparse matrix assembly  and vector assembly
export AbstractSysmatAssembler, SysmatAssemblerSparse, SysmatAssemblerSparseSymm, SysmatAssemblerSparseHRZLumpingSymm, startassembly!, assemble!, makematrix!, AbstractSysvecAssembler, SysvecAssembler, startassembly!, assemble!, makevector!

using .IntegRuleModule: AbstractIntegRule, TriRule, GaussRule, TetRule, PointRule, SimplexRule, TrapezoidalRule, NodalSimplexRule, NodalTensorProductRule
# Exported: type for various integration rules
export AbstractIntegRule, TriRule, GaussRule, TetRule, PointRule, SimplexRule, TrapezoidalRule, NodalSimplexRule, NodalTensorProductRule

using .IntegDomainModule: IntegDomain, otherdimensionunity, Jacobianpoint, Jacobiancurve, Jacobiansurface, Jacobianvolume, Jacobianmdim, integrationdata
# Exported: type to handle  integration data for various manifold dimensions
export IntegDomain, otherdimensionunity, Jacobianpoint, Jacobiancurve, Jacobiansurface, Jacobianvolume, Jacobianmdim, integrationdata

using .FEMMBaseModule: AbstractFEMM, FEMMBase, associategeometry!, integratefieldfunction, integratefunction, transferfield!, distribloads, connectionmatrix, fieldfromintegpoints, elemfieldfromintegpoints, innerproduct
# Exported: type base discretization methods
export AbstractFEMM, FEMMBase, associategeometry!, integratefieldfunction, integratefunction, transferfield!, distribloads, connectionmatrix, fieldfromintegpoints, elemfieldfromintegpoints, innerproduct

###########################################################################
# Mesh-generation functionality for various shapes
###########################################################################
using .MeshQuadrilateralModule: Q4annulus, Q4quadrilateral, Q4elliphole, Q4block, Q4blockx, Q4refine, Q8block, Q4toQ8, Q8annulus, Q8blockx, Q4spheren, Q4circlen
# Exported: mesh generation functions for quadrilateral elements
export Q4annulus, Q4quadrilateral, Q4elliphole, Q4block, Q4blockx, Q4refine, Q8block, Q4toQ8, Q8annulus, Q8blockx, Q4spheren, Q4circlen

using .MeshLineModule: L2block, L2blockx, L3blockx
# Exported: mesh generation functions for line elements
export L2block, L2blockx, L3blockx

using .MeshTriangleModule: T3blockx, T3block,  T3toT6,  T6block,  Q4toT3,  T3refine, T6blockx, T3annulus, T6annulus
# Exported: mesh generation functions for triangular elements
export  T3blockx, T3block,  T3toT6,  T6block,  Q4toT3,  T3refine, T6blockx, T3annulus, T6annulus

using .MeshHexahedronModule: H8block,  H8blockx,  H8sphere,  H8refine, H8hexahedron, H8extrudeQ4, H8spheren, H8voximg,  H8layeredplatex, H8elliphole, H8toH27,  H27block, H20block,  H8toH20, H20blockx, H27blockx, H8cylindern, T4toH8
# Exported: mesh generation functions for hexahedral elements
export  H8block,  H8blockx,  H8sphere,  H8refine, H8hexahedron, H8extrudeQ4, H8spheren, H8voximg,  H8layeredplatex, H8elliphole, H8toH27,  H27block, H20block,  H8toH20, H20blockx, H27blockx, H8cylindern, T4toH8

using .MeshTetrahedronModule: T4block, T4blockx, T4toT10, T10toT4, T10block, T10blockx, T10layeredplatex, T4meshedges, T4voximg, T4refine, T10refine, T4refine20, T4quartercyln
# Exported: mesh generation functions for tetrahedral elements
export  T4block, T4blockx, T4toT10, T10toT4, T10block, T10blockx, T10layeredplatex, T4meshedges, T4voximg, T4refine, T10refine, T4refine20, T4quartercyln

###########################################################################
# Abstract material
###########################################################################
using .MatModule: AbstractMat, massdensity
# Exported: abstract type of material
export AbstractMat, massdensity

end
