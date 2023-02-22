"""
    MeshMixinsModule

Module  for generation of various meshes.
"""
module MeshMixinsModule

__precompile__(true)

import LinearAlgebra: norm
using ..FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import ..FESetModule: AbstractFESet, FESetT3, FESetT6, FESetQ4, connasarray, subset
import ..FENodeSetModule: FENodeSet
import ..MeshUtilModule: makecontainer, addhyperface!, findhyperface!, linearspace
import ..MeshModificationModule:
    meshboundary,
    connectednodes,
    fusenodes,
    updateconn!,
    compactnodes,
    renumberconn!,
    meshsmoothing
import ..MeshSelectionModule: selectelem, findunconnnodes
import ..MeshQuadrilateralModule: Q4block
import ..MeshTriangleModule: Q4toT3, T3toT6
import Statistics: mean

"""
    T3blockrand(Length::FFlt, Width::FFlt, nL::FInt, nW::FInt)

T3 mesh of a rectangle with random orientation.
"""
function T3blockrand(Length::FFlt, Width::FFlt, nL::FInt, nW::FInt)
    # nL = Int(round(nL/2))
    # nW = Int(round(nW/2))
    fens, fes = Q4block(Length, Width, nL, nW) # Mesh
    fens, fes = Q4toT3(fens, fes, :random) # Mesh
    # fens,fes = T3toQ4(fens,fes); # Mesh
    # fens,fes = Q4toT3(fens,fes); # Mesh
    # File = "mesh_a.vtk"
    # VTK.vtkexportmesh(File, fens, fes)

    fixedv = fill(0, count(fens))
    fixedv[connectednodes(meshboundary(fes))] .= 1
    fens = meshsmoothing(fens, fes; fixedv = fixedv)
    # File = "mesh_b.vtk"
    # VTK.vtkexportmesh(File, fens, fes)
    return fens, fes
end

"""
    T3toQ4(fens,fes)

Convert a T3 triangular mesh into a quadrilateral mesh.

Convert a mesh of triangles T3 to three quadrilateral Q4 each.
"""
function T3toQ4(fens, fes)
    fens, fes = T3toT6(fens, fes)
    onfens = count(fens)
    xyz1 = deepcopy(fens.xyz)
    id = zeros(count(fens) + count(fes))
    id[1:count(fens)] = 1:count(fens)
    xyz = zeros(count(fens) + count(fes), 3)
    xyz[1:count(fens), 1:size(xyz1, 2)] = xyz1
    # construct new geometry cells
    nconns = zeros(FInt, 3 * count(fes), 4)
    conns = connasarray(fes)
    nc = onfens + 1
    gc = 1
    for i in 1:count(fes)
        conn = conns[i, :]
        id[nc] = nc
        xyz[nc, :] = mean(xyz[conn[1:3], :], dims = 1)
        nconns[gc, :] = [conn[1] conn[4] nc conn[6]]
        gc = gc + 1
        nconns[gc, :] = [conn[2] conn[5] nc conn[4]]
        gc = gc + 1
        nconns[gc, :] = [conn[3] conn[6] nc conn[5]]
        gc = gc + 1
        nc = nc + 1
    end
    fes = FESetQ4(nconns)
    fens.xyz = xyz
    return fens, fes
end


end
