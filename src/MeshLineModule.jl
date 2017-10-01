"""
    MeshLineModule

Module  for generation of meshes composed of line (curve) elements.
"""
module MeshLineModule

export L2block, L2blockx, L3blockx

using FinEtools.FTypesModule
using FinEtools.FESetModule
using FinEtools.FENodeSetModule
using FinEtools.MeshModificationModule

"""
L2block(Length::FFlt, nL::FInt)

Mesh of a 1-D block of L2 finite elements.
"""
function L2block(Length::FFlt, nL::FInt)
    fens,fes = L2blockx(collect(squeeze(linspace(0,Length,nL+1)',1)));
end


"""
L2blockx(xs::FFltVec)

Graded mesh of a 1-D block, L2 finite elements.
"""
function L2blockx(xs::FFltVec)
    xyz =reshape(sort(xs),length(xs),1);
    ncells=length(xs)-1;

    # create the nodes
    fens = FENodeSetModule.FENodeSet(xyz);
    # Create the finite elements
    fes = FESetModule.FESetL2([(1:ncells) (2:ncells+1)]);

    return fens,fes;
end


"""
    L3blockx(xs::FFltVec)

Graded mesh of a 1-D block, L2 finite elements.
"""
function L3blockx(xs::FFltVec)
    fens,fes = L2blockx(xs)
    nxyz = zeros(count(fes), size(fens.xyz, 2))
    nconn = zeros(FInt, count(fes), 3)
    N = count(fens)
    for i = 1:count(fes)
        N = N + 1
        nxyz[i, :] = mean(fens.xyz[fes.conn[i, :], :], 1)
        nconn[i, :] = vcat(fes.conn[i, :], [N])
    end
    fens = FENodeSetModule.FENodeSet([fens.xyz; nxyz]);
    # Create the finite elements
    fes = FESetModule.FESetL3(nconn);

    return fens,fes;
end


end
