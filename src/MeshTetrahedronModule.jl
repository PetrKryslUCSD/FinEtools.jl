"""
    MeshTetrahedronModule

Module  for generation of meshes composed of tetrahedra.
"""
module MeshTetrahedronModule

export  T4block, T4blockx, T4toT10, T10block, T10blockx, T10compositeplatex

using FinEtools.FTypesModule
using FinEtools.FESetModule
using FinEtools.FENodeSetModule
using FinEtools.MeshUtilModule
using FinEtools.MeshModificationModule
using FinEtools.MeshSelectionModule

"""
    T4block(Length::FFlt, Width::FFlt, Height::FFlt,
       nL::FInt, nW::FInt, nH::FInt, orientation::Symbol)

Generate a tetrahedral mesh  of the 3D block.

Four-node tetrahedra in a regular arrangement, with uniform spacing between
the nodes, with a given orientation of the diagonals.

The mesh is produced by splitting each logical  rectangular cell into six
tetrahedra.
Range =<0, Length> x <0, Width> x <0, Height>
Divided into elements: nL,  nW,  nH in the first,  second,  and
third direction (x, y, z).
"""
function T4block(Length::FFlt, Width::FFlt, Height::FFlt,
  nL::FInt, nW::FInt, nH::FInt, orientation::Symbol=:a)
    return T4blockx(collect(linspace(0.0, Length, nL+1)),
                    collect(linspace(0.0, Width, nW+1)),
                    collect(linspace(0.0, Height, nH+1)), orientation);
end

"""
    T4blockx(xs::FFltMat, ys::FFltMat, zs::FFltMat, orientation::Symbol)

Generate a graded tetrahedral mesh  of a 3D block.

Four-node tetrahedra in a regular arrangement, with non-uniform given spacing
between the nodes, with a given orientation of the diagonals.

The mesh is produced by splitting each logical  rectangular cell into six
tetrahedra.
"""
function T4blockx(xs::FFltMat, ys::FFltMat, zs::FFltMat, orientation::Symbol)
    return T4blockx(vec(xs), vec(ys), vec(zs), orientation)
end

"""
    T4blockx(xs::FFltVec, ys::FFltVec, zs::FFltVec, orientation::Symbol)

Generate a graded tetrahedral mesh  of a 3D block.

Four-node tetrahedra in a regular arrangement, with non-uniform given spacing
between the nodes, with a given orientation of the diagonals.

The mesh is produced by splitting each logical  rectangular cell into six
tetrahedra.
"""
function T4blockx(xs::FFltVec, ys::FFltVec, zs::FFltVec, orientation::Symbol)
    nL =length(xs)-1;
    nW =length(ys)-1;
    nH =length(zs)-1;
    nnodes=(nL+1)*(nW+1)*(nH+1);
    ncells=6*(nL)*(nW)*(nH);
    xyzs=zeros(FFlt, nnodes, 3);
    conns=zeros(FInt, ncells, 4);
    if (orientation==:a)
        t4ia = [1  8  5  6; 3  4  2  7; 7  2  6  8; 4  7  8  2; 2  1  6  8; 4  8  1  2];
        t4ib = [1  8  5  6; 3  4  2  7; 7  2  6  8; 4  7  8  2; 2  1  6  8; 4  8  1  2];
    elseif (orientation==:b)
        t4ia = [2 7 5 6; 1 8 5 7; 1 3 4 8; 2 1 5 7; 1 2 3 7; 3 7 8 1];
        t4ib = [2 7 5 6; 1 8 5 7; 1 3 4 8; 2 1 5 7; 1 2 3 7; 3 7 8 1];
    elseif (orientation==:ca)
        t4ia = [8  4  7  5; 6  7  2  5; 3  4  2  7; 1  2  4  5; 7  4  2  5];
        t4ib = [7  3  6  8; 5  8  6  1; 2  3  1  6; 4  1  3  8; 6  3  1  8];
    elseif (orientation==:cb)
        t4ia = [7  3  6  8; 5  8  6  1; 2  3  1  6; 4  1  3  8; 6  3  1  8];
        t4ib = [8  4  7  5; 6  7  2  5; 3  4  2  7; 1  2  4  5; 7  4  2  5];
      else
        error("Unknown orientation")
      end
    f=1;
    for k=1:(nH+1)
        for j=1:(nW+1)
            for i=1:(nL+1)
                xyzs[f, :]=[xs[i] ys[j] zs[k]];
                f=f+1;
            end
        end
    end

    fens=FENodeSetModule.FENodeSet(xyzs);

    function node_numbers(i::FInt, j::FInt, k::FInt, nL::FInt, nW::FInt, nH::FInt)
        f=(k-1)*((nL+1)*(nW+1))+(j-1)*(nL+1)+i;
        nn=[f (f+1)  f+(nL+1)+1 f+(nL+1)];
        return [nn broadcast(+, nn, (nL+1)*(nW+1))];
    end

    gc=1;
    for i=1:nL
        for j=1:nW
            for k=1:nH
                nn=node_numbers(i, j, k, nL, nW, nH);
                if (mod(sum( [i, j, k] ), 2)==0)
                    t4i =t4ib;
                else
                    t4i =t4ia;
                end
                for r=1:size(t4i, 1)
                    for c1=1:size(t4i, 2)
                        conns[gc, c1] = nn[t4i[r, c1]];
                    end
                    gc=gc+1;
                end
            end
        end
    end
    fes = FESetModule.FESetT4(conns[1:gc-1, :]);

    return fens, fes
end

"""
    T4toT10(fens::FENodeSetModule.FENodeSet,  fes::FESetModule.FESetT4)

Convert a mesh of tetrahedra of type T4 (four-node) to tetrahedra T10.
"""
function  T4toT10(fens::FENodeSetModule.FENodeSet,  fes::FESetModule.FESetT4)
    nedges=6;
    ec = [1  2; 2  3; 3  1; 4  1; 4  2; 4  3];
    # Additional node numbers are numbered from here
    newn = FENodeSetModule.count(fens)+1;
    # make a search structure for edges
    edges = MeshUtilModule.makecontainer();
    for i= 1:size(fes.conn, 1)
        conn = fes.conn[i, :];
        for J = 1:nedges
            ev=conn[ec[J, :]];
            newn = MeshUtilModule.addhyperface!(edges,  ev,  newn);
        end
    end
    xyz1 =fens.xyz;             # Pre-existing nodes
    # Allocate for vertex nodes plus edge nodes plus face nodes
    xyz =zeros(FFlt, newn-1, 3);
    xyz[1:size(xyz1, 1), :] = xyz1; # existing nodes are copied over
    # calculate the locations of the new nodes
    # and construct the new nodes
    for i in keys(edges)
        C=edges[i];
        for J = 1:length(C)
          ix = vec([item for item in C[J].o])
          push!(ix,  i) # Add the anchor point as well
          xyz[C[J].n, :] = mean(xyz[ix, :], 1);
        end
    end
    fens =FENodeSetModule.FENodeSet(xyz);
    # construct new geometry cells
    nconn=zeros(FInt, size(fes.conn, 1), 10);
    nc=1;
    for i= 1:size(fes.conn, 1)
        conn = fes.conn[i, :];
        econn=zeros(FInt, 1, nedges);
        for J = 1:nedges
            ev=conn[ec[J, :]];
            h, n=MeshUtilModule.findhyperface!(edges,  ev);
            econn[J]=n;
        end
        nconn[nc, :] = vcat(vec(conn),   vec(econn))
        nc= nc+ 1;
    end
    fes = FESetModule.FESetT10(nconn);
    return fens, fes;
end

"""
    T10block(Length::FFlt, Width::FFlt, Height::FFlt,
      nL::FInt, nW::FInt, nH::FInt; orientation::Symbol=:a)

Generate a tetrahedral  mesh of T10 elements  of a rectangular block.
"""
function T10block(Length::FFlt, Width::FFlt, Height::FFlt,
    nL::FInt, nW::FInt, nH::FInt; orientation::Symbol=:a)
    fens, fes = T4block(Length, Width, Height, nL, nW, nH, orientation);
    fens, fes = T4toT10(fens, fes);
    return fens, fes
end

"""
    T10blockx(xs::FFltMat, ys::FFltMat, zs::FFltMat, orientation::Symbol = :a)

Generate a graded 10-node tetrahedral mesh  of a 3D block.

10-node tetrahedra in a regular arrangement, with non-uniform given spacing
between the nodes, with a given orientation of the diagonals.

The mesh is produced by splitting each logical  rectangular cell into six
tetrahedra.
"""
function T10blockx(xs::FFltMat, ys::FFltMat, zs::FFltMat, orientation::Symbol = :a)
    fens, fes =  T4blockx(vec(xs), vec(ys), vec(zs), orientation)
    fens, fes = T4toT10(fens, fes);
    return fens, fes
end

function T10blockx(xs::FFltVec, ys::FFltVec, zs::FFltVec, orientation::Symbol = :a)
    fens, fes =  T4blockx(vec(xs), vec(ys), vec(zs), orientation)
    fens, fes = T4toT10(fens, fes);
    return fens, fes
end

"""
    T10compositeplatex(xs::FFltVec, ys::FFltVec, ts::FFltVec, nts::FIntVec,
        orientation::Symbol = :a)

T10 mesh for a layered block (composite plate) with specified in plane coordinates.

xs,ys =Locations of the individual planes of nodes.
ts= Array of layer thicknesses,
nts= array of numbers of elements per layer

The finite elements of each layer are labeled with the layer number, starting
from 1.
"""
function T10compositeplatex(xs::FFltVec, ys::FFltVec, ts::FFltVec, nts::FIntVec,
    orientation::Symbol = :a)
    tolerance = minimum(abs.(ts))/maximum(nts)/10.;
    @assert length(ts) >= 1
    @assert sum(nts) >= length(ts)
    zs = collect(linspace(0.0, ts[1], nts[1]+1))
    for layer = 2:length(ts)
        oz = collect(linspace(sum(ts[1:layer-1]), sum(ts[1:layer]), nts[layer]+1))
        zs = vcat(zs, oz[2:end])
    end
    fens, fes = T10blockx(xs, ys, zs, orientation);
    List = selectelem(fens, fes, box = [-Inf Inf -Inf Inf 0.0 ts[1]],
            inflate = tolerance)
    fes.label[List] = 1
    for layer = 2:length(ts)
        List = selectelem(fens, fes,
            box = [-Inf Inf -Inf Inf sum(ts[1:layer-1]) sum(ts[1:layer])],
            inflate = tolerance)
        fes.label[List] = layer
    end
    return fens,fes
end

end
