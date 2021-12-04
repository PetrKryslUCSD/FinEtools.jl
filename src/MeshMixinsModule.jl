"""
    MeshMixinsModule

Module  for generation of various meshes.
"""
module MeshMixinsModule

__precompile__(true)

import LinearAlgebra: norm
using ..FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import ..FESetModule: AbstractFESet, FESetT3, FESetT6, FESetQ4, connasarray, subset
import ..FENodeSetModule: FENodeSet
import ..MeshUtilModule: makecontainer, addhyperface!, findhyperface!, linearspace
import ..MeshModificationModule: meshboundary, connectednodes, fusenodes, updateconn!, compactnodes, renumberconn!, meshsmoothing
import ..MeshSelectionModule: selectelem, findunconnnodes
import ..MeshQuadrilateralModule: Q4block
import ..MeshTriangleModule: Q4toT3
import Statistics: mean

"""
    T3blockrand(Length::FFlt, Width::FFlt, nL::FInt, nW::FInt)

T3 mesh of a rectangle with random orientation.
"""
function T3blockrand(Length::FFlt, Width::FFlt, nL::FInt, nW::FInt)
    fens,fes = Q4block(Length,Width,nL,nW); # Mesh
    fens,fes = Q4toT3(fens,fes,:random); # Mesh
    # File = "mesh_a.vtk"
    # VTK.vtkexportmesh(File, fens, fes)

    fixedv = fill(0, count(fens))
    fixedv[connectednodes(meshboundary(fes))] .= 1
    fens = meshsmoothing(fens, fes; fixedv = fixedv)
    # File = "mesh_b.vtk"
    # VTK.vtkexportmesh(File, fens, fes)
    return fens,fes
end



end
