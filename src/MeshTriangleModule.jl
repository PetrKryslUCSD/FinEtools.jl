"""
    MeshTriangleModule

Module  for generation of meshes composed of triangles.
"""
module MeshTriangleModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FESetModule: AbstractFESet, FESetT3, FESetT6, FESetQ4, connasarray
import FinEtools.FENodeSetModule: FENodeSet
import FinEtools.MeshUtilModule: makecontainer, addhyperface!, findhyperface!, linearspace
import Statistics: mean

"""
    T3blockx(xs::FFltVec, ys::FFltVec, orientation::Symbol=:a)

T3 Mesh of a rectangle.
"""
function T3blockx(xs::FFltVec, ys::FFltVec, orientation::Symbol=:a)
    if (orientation!=:a) && (orientation!=:b)
        error("Cannot handle orientation : $orientation")
    end
    nL =length(xs)-1;
    nW =length(ys)-1;
    nnodes=(nL+1)*(nW+1);
    ncells=2*(nL)*(nW);
    xys=zeros(FFlt,nnodes,2);
    conns=zeros(FInt,ncells,3);
    f=1;
    for j=1:(nW+1)
        for i=1:(nL+1)
            xys[f, 1] = xs[i]
            xys[f, 2] = ys[j]
            f=f+1;
        end
    end
    fens = FENodeSet(xys);

    gc=1;
    for i=1:nL
        for j=1:nW
            f=(j-1)*(nL+1)+i;
            if     (orientation==:a)
                nn=[f (f+1) f+(nL+1)];
            elseif (orientation==:b)
                nn=[f (f+1) f+(nL+1)+1];
            end
            conns[gc,:]=nn;
            gc=gc+1;
            if     (orientation==:a)
                nn=[(f+1)  f+(nL+1)+1 f+(nL+1)];
            elseif (orientation==:b)
                nn=[f  f+(nL+1)+1 f+(nL+1)];
            end
            conns[gc,:]=nn;
            gc=gc+1;
        end
    end
    fes = FESetT3(conns);
    return fens,fes
end

function T3blockx(xs::FFltMat,ys::FFltMat,orientation::Symbol=:a)
    return T3blockx(vec(xs),vec(ys),orientation)
end

"""
    T3block(Length::FFlt, Width::FFlt, nL::FInt, nW::FInt, orientation::Symbol=:a)

T3 Mesh of a rectangle.
"""
function T3block(Length::FFlt, Width::FFlt, nL::FInt, nW::FInt, orientation::Symbol=:a)
    return T3blockx(FFltVec(linearspace(0.0,Length,nL+1)), FFltVec(linearspace(0.0,Width,nW+1)), orientation)
end

"""
    T3toT6(fens::FENodeSet, fes::FESetT3)

Convert a mesh of triangle T3 (three-node) to triangle T6.
"""
function T3toT6(fens::FENodeSet, fes::FESetT3)
    nedges=3;
    ec = [1 2; 2 3; 3 1];
    conns = connasarray(fes);
    # Additional node numbers are numbered from here
    newn=count(fens)+1;
    # make a search structure for edges
    edges=makecontainer();
    for i= 1:size(conns,1)
        conn = conns[i,:];
        for J = 1:nedges
            ev=conn[ec[J,:]];
            newn = addhyperface!(edges, ev, newn);
        end
    end
    xyz1 =fens.xyz;             # Pre-existing nodes
    # Allocate for vertex nodes plus edge nodes plus face nodes
    xyz =zeros(FFlt,newn-1,size(fens.xyz,2));
    xyz[1:size(xyz1,1),:] = xyz1; # existing nodes are copied over
    # calculate the locations of the new nodes
    # and construct the new nodes
    for i in keys(edges)
        C=edges[i];
        for J = 1:length(C)
            ix = vec([item for item in C[J].o])
            push!(ix,  i)
            xyz[C[J].n, :] = mean(xyz[ix, :], dims = 1);
        end
    end
    # construct new geometry cells
    nconns =zeros(FInt,size(conns,1),6);
    nc=1;
    for i= 1:size(conns,1)
        conn = conns[i,:];
        econn=zeros(FInt,1,nedges);
        for J = 1:nedges
            ev=conn[ec[J,:]];
            h,n=findhyperface!(edges, ev);
            econn[J]=n;
        end
        nconns[nc,:] =vcat(vec(conn), vec(econn));
        nc= nc+ 1;
    end
    fens =FENodeSet(xyz);
    fes = FESetT6(nconns) ;
    return fens,fes;
end

"""
    T6block(Length::FFlt, Width::FFlt, nL::FInt, nW::FInt, orientation::Symbol=:a)

Mesh of a rectangle of T6 elements.
"""
function T6block(Length::FFlt, Width::FFlt, nL::FInt, nW::FInt, orientation::Symbol=:a)
    fens,fes = T3block(Length,Width,nL,nW,orientation);
    fens,fes = T3toT6(fens,fes);
end

"""
    T6blockx(xs::FFltVec, ys::FFltVec, orientation::Symbol=:a)

Graded mesh of a 2-D block of T6 finite elements.
"""
function T6blockx(xs::FFltVec, ys::FFltVec, orientation::Symbol=:a)
    fens, fes = T3blockx(xs, ys, orientation);
    fens, fes = T3toT6(fens, fes);
end

"""
    Q4toT3(fens::FENodeSet, fes::FESetQ4, orientation::Symbol=:default)

Convert a mesh of quadrilateral Q4's to two T3 triangles  each.
"""
function Q4toT3(fens::FENodeSet, fes::FESetQ4, orientation::Symbol=:default)
    connl1=[1  2  3];
    connl2=[1  3  4];
    if orientation==:alternate
        connl1=[1, 2, 4];
        connl2=[3, 4, 2];
    end
    nedges=4;
    nconns=zeros(FInt,2*count(fes),3);
    conns = connasarray(fes)
    nc=1;
    for i= 1:size(conns, 1)
        conn = conns[i,:];
        nconns[nc,:] =conn[connl1];
        nc= nc+ 1;
        nconns[nc,:] =conn[connl2];
        nc= nc+ 1;
    end
    nfes = FESetT3(nconns);
    return fens,nfes            # I think I should not be overwriting the input!
end

"""
    T3refine(fens::FENodeSet,fes::FESetT3)

Refine a mesh of 3-node triangles by quadrisection.
"""
function T3refine(fens::FENodeSet, fes::FESetT3)
    fens,fes = T3toT6(fens,fes);
    nconn=zeros(FInt,4*size(fes.conn,1),3);
    conns = connasarray(fes)
    nc=1;
    for i= 1:size(conns,1)
        c=conns[i,:];
        nconn[nc,:] =c[[1,4,6]];        nc= nc+ 1;
        nconn[nc,:] =c[[2,5,4]];        nc= nc+ 1;
        nconn[nc,:] =c[[3,6,5]];        nc= nc+ 1;
        nconn[nc,:] =c[[4,5,6]];        nc= nc+ 1;
    end
    nfes = FESetT3(nconn);
    return fens,nfes            # I think I should not be overwriting the input!
end

"""
    T3annulus(rin::FFlt, rex::FFlt, nr::FInt, nc::FInt, Angl::FFlt, orientation::Symbol=:a)

Mesh of an annulus segment.

Mesh of an annulus segment, centered at the origin, with internal radius `rin`,
and  external radius `rex`, and  development angle `Angl` (in radians). Divided
into elements: nr, nc in the radial and circumferential direction respectively.
"""
function T3annulus(rin::FFlt, rex::FFlt, nr::FInt, nc::FInt, Angl::FFlt, orientation::Symbol=:a)
    trin=min(rin,rex);
    trex=max(rin,rex);
    fens,fes = T3block(trex-trin,Angl,nr,nc,orientation);
    xy=fens.xyz;
    for i=1:count(fens)
        r=trin+xy[i,1]; a=xy[i,2];
        xy[i,:]=[r*cos(a) r*sin(a)];
    end
    fens.xyz=xy;
    return fens,fes
end

"""
    T6annulus(rin::FFlt, rex::FFlt, nr::FInt, nc::FInt, Angl::FFlt, orientation::Symbol=:a)

Mesh of an annulus segment.

Mesh of an annulus segment, centered at the origin, with internal radius `rin`,
and  external radius `rex`, and  development angle `Angl` (in radians). Divided
into elements: nr, nc in the radial and circumferential direction respectively.
"""
function T6annulus(rin::FFlt, rex::FFlt, nr::FInt, nc::FInt, Angl::FFlt, orientation::Symbol=:a)
    trin=min(rin,rex);
    trex=max(rin,rex);
    fens,fes = T3block(trex-trin,Angl,nr,nc,orientation);
    fens, fes = T3toT6(fens, fes)
    xy=fens.xyz;
    for i=1:count(fens)
        r=trin+xy[i,1]; a=xy[i,2];
        xy[i,:]=[r*cos(a) r*sin(a)];
    end
    fens.xyz=xy;
    return fens,fes
end

end
