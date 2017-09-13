"""
    MeshQuadrilateralModule

Module  for generation of meshes composed of quadrilaterals.
"""
module MeshQuadrilateralModule

export  Q4annulus,  Q4quadrilateral,  Q4elliphole,  Q4block,  Q4blockx,  Q4refine
export Q8block,  Q4toQ8,  Q8annulus

using FinEtools.FTypesModule
using FinEtools.FESetModule
using FinEtools.FENodeSetModule
using FinEtools.MeshModificationModule
using FinEtools.MeshUtilModule

"""
    Q4annulus(rin::FFlt, rex::FFlt, nr::FInt, nc::FInt, Angl::FFlt)

Mesh of an annulus segment.

Mesh of an annulus segment, centered at the origin, with internal radius `rin`,
and  external radius `rex`, and  development angle `Angl` (in radians). Divided
into elements: nr, nc in the radial and circumferential direction respectively.
"""
function Q4annulus(rin::FFlt, rex::FFlt, nr::FInt, nc::FInt, Angl::FFlt)
    trin=min(rin,rex);
    trex=max(rin,rex);
    fens,fes =Q4block(trex-trin,Angl,nr,nc);
    xy=fens.xyz;
    for i=1:FENodeSetModule.count(fens)
        r=trin+xy[i,1]; a=xy[i,2];
        xy[i,:]=[r*cos(a) r*sin(a)];
    end
    fens.xyz=xy;
    return fens,fes
end

"""
    Q4quadrilateral(xyz::FFltMat, nL::FInt, nW::FInt)

Mesh of a general quadrilateral given by the location of the vertices.
"""
function Q4quadrilateral(xyz::FFltMat, nL::FInt, nW::FInt)
    npts=size(xyz,1);
    if npts==2
        lo=minimum(xyz,1);
        hi=maximum(xyz,1);
        xyz=[[lo[1] lo[2]];
        [hi[1] lo[2]];
        [hi[1] hi[2]];
        [lo[1] hi[2]]];
    elseif npts!=4
        error("Need 2 or 4 points");
    end

    fens,fes = Q4block(2.,2.,nL,nW);

    xyz1=fens.xyz;
    if (size(xyz1,2)<size(xyz,2))
        nxyz1=zeros(FFlt,size(xyz1,1),size(xyz,2));
        nxyz1[:,1:size(xyz1,2)]=xyz1;
        xyz1=nxyz1;
    end

    dummy = FESetModule.FESetQ4(reshape(collect(1:4),1,4))
    pxyz = fens.xyz;
    for i = 1:FENodeSetModule.count(fens)
        N = FESetModule.bfun(dummy, broadcast(-, pxyz[i,:], 1.0));# shift coordinates by -1
        pxyz[i,:] = N'*xyz;
    end
    copy!(fens.xyz, xyz1);
    return fens,fes
end

"""
    Q4elliphole(xradius::FFlt, yradius::FFlt, L::FFlt, H::FFlt,
      nL::FInt, nH::FInt, nW::FInt)

Mesh of one quarter of a rectangular plate with an elliptical hole.

    xradius,yradius = radius of the ellipse,
    L,H= and dimensions of the plate,
    nL,nH= numbers of edges along the side of the plate; this also happens
      to be the number of edges along the circumference of the elliptical
      hole
    nW= number of edges along the remaining straight edge (from the hole
      in the direction of the length),
"""
function Q4elliphole(xradius::FFlt, yradius::FFlt, L::FFlt, H::FFlt,
    nL::FInt, nH::FInt, nW::FInt)
    dA =pi/2/(nL +nH);
    tolerance = (xradius+yradius)/(nL*nH)/100;
    fens= nothing; fes= nothing;
    for i= 1:nH
        xy = [xradius*cos((i-1)*dA) yradius*sin((i-1)*dA);
        L (i-1)/nH*H;
        L (i)/nH*H;
        xradius*cos((i)*dA) yradius*sin((i)*dA)];
        fens1,fes1 = Q4quadrilateral(xy,nW,1);
        if (fens == nothing)
            fens = fens1; fes = fes1;
        else
            fens,fes1,fes2 = MeshModificationModule.mergemeshes(fens1, fes1, fens, fes, tolerance);
            fes = FESetModule.cat(fes1,fes2);
        end
    end
    for i= 1:nL
        xy = [xradius*cos((nH+i-1)*dA)   yradius*sin((nH+i-1)*dA);
        (nL-i+1)/nL*L   H;
        (nL-i)/nL*L  H;
        xradius*cos((nH+i)*dA)   yradius*sin((nH+i)*dA)];
        fens1,fes1 = Q4quadrilateral(xy,nW,1);
        fens,fes1,fes2 = MeshModificationModule.mergemeshes(fens1, fes1, fens, fes, tolerance);
        fes = FESetModule.cat(fes1, fes2);
    end
    return fens,fes
end

"""
    Q4block(Length::FFlt, Width::FFlt, nL::FInt, nW::FInt)

Mesh of a rectangle, Q4 elements.

Divided into elements: nL, nW in the first, second (x,y).
"""
function Q4block(Length::FFlt, Width::FFlt, nL::FInt, nW::FInt)
    return Q4blockx(collect(linspace(0.0,Length,nL+1)),
        collect(linspace(0.0,Width,nW+1)));
end

"""
    Q4blockx(xs::FFltVec, ys::FFltVec)

Graded mesh  of a rectangle, Q4 finite elements.

    Mesh of a 2-D block, Q4 finite elements. The nodes are located at the
    Cartesian product of the two intervals on the input.  This allows for
    construction of graded meshes.

    xs,ys - Locations of the individual planes of nodes.
"""
function Q4blockx(xs::FFltVec, ys::FFltVec)
    nL = length(xs) - 1;
    nW = length(ys) - 1;

    nnodes = (nL+1) * (nW+1);
    ncells = nL * nW;

    # preallocate node locations
    xyz = zeros(FFlt, nnodes, 2);
    k = 1;
    for j = 1:(nW+1)
        for i = 1:(nL+1)
            xyz[k,:] = [xs[i] ys[j]];
            k = k + 1;
        end
    end
    # create the nodes
    fens = FENodeSetModule.FENodeSet(xyz);

    #preallocate connectivity matrix
    conn = zeros(FInt, ncells, 4);

    function  nodenumbers(i,j,nL,nW)
        f = (j-1) * (nL+1) + i;
        nn = [f, (f+1), f+(nL+1)+1, f+(nL+1)];
        return nn
    end

    k = 1;
    for i = 1:nL
        for j = 1:nW
            conn[k,:] = nodenumbers(i,j,nL,nW);
            k = k + 1;
        end
    end
    # create the cells
    fes = FESetModule.FESetQ4(conn);

    return fens,fes;
end

function Q4blockx(xs::AbstractVector, ys::AbstractVector)
    return Q4blockx(FFltVec(xs), FFltVec(ys))
end

"""
    Q8block(Length::FFlt, Width::FFlt, nL::FInt, nW::FInt)

Mesh of a rectangle of Q8 elements.
"""
function Q8block(Length::FFlt, Width::FFlt, nL::FInt, nW::FInt)
    fens,fes  = Q4block(Length, Width, nL, nW);
    fens,fes = Q4toQ8(fens, fes);
end

"""
    Q4toQ8(fens::FENodeSet, fes::FESetQ4)

Convert a mesh of quadrilateral Q4 to quadrilateral Q8.
"""
function Q4toQ8(fens::FENodeSet, fes::FESetQ4)
    nedges=4;
    ec = [1  2; 2  3; 3  4; 4  1];
    conns = fes.conn;
    # Additional node numbers are numbered from here
    newn = count(fens)+1;
    # make a search structure for edges
    edges = MeshUtilModule.makecontainer();
    for i= 1:size(conns,1)
        conn = conns[i,:];
        for J = 1:nedges
            ev = conn[ec[J,:]];
            newn = MeshUtilModule.addhyperface!(edges, ev, newn);
        end
    end
    xyz1 =fens.xyz;             # Pre-existing nodes
    # Allocate for vertex nodes plus edge nodes plus face nodes
    xyz =zeros(FFlt,newn-1,size(xyz1,2));
    xyz[1:size(xyz1,1),:] = xyz1; # existing nodes are copied over
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
    # construct new geometry cells
    nconns =zeros(FInt,size(conns,1),8);
    nc=1;
    for i= 1:size(conns,1)
        conn = conns[i,:];
        econn=zeros(FInt,1,nedges);
        for J = 1:nedges
            ev=conn[ec[J,:]];
            h,n = MeshUtilModule.findhyperface!(edges, ev);
            econn[J]=n;
        end
        nconns[nc,:] =vcat(vec(conn), vec(econn));
        nc= nc+ 1;
    end
    fens = FENodeSet(xyz);
    fes = FESetQ8(nconns);
    return fens, fes
end

"""
    Q4refine(fens::FENodeSet, fes::FESetQ4)

Refine a mesh of quadrilaterals by bisection.
"""
function Q4refine(fens::FENodeSet, fes::FESetQ4)
    nedges=4;
    ec = [1  2; 2  3; 3  4; 4  1];
    # make a search structure for edges
    # Additional node numbers are numbered from here
    newn=FENodeSetModule.count(fens)+1;
    # make a search structure for edges
    edges=MeshUtilModule.makecontainer();
    for i= 1:size(fes.conn,1)
        conn = fes.conn[i,:];
        for J = 1:nedges
            ev=conn[ec[J,:]];
            newn = MeshUtilModule.addhyperface!(edges, ev, newn);
        end
    end
    newn=  newn+size(fes.conn,1) # add the interior nodes to the total
    xyz1 =fens.xyz;             # Pre-existing nodes
    # Allocate for vertex nodes plus edge nodes plus face nodes
    xyz =zeros(FFlt,newn-1,size(xyz1,2));
    xyz[1:size(xyz1,1),:] = xyz1; # existing nodes are copied over
    # calculate the locations of the new nodes
    # and construct the new nodes
    for i in keys(edges)
        C=edges[i];
        for J = 1:length(C)
            ix = vec([item for item in C[J].o])
            push!(ix, i)
            xyz[C[J].n,:] = mean(xyz[ix,:],1);
        end
    end
    # construct new geometry cells: for new elements out of one old one
    nconn =zeros(FInt,4*size(fes.conn,1),4);
    nc=1;
    for i= 1:size(fes.conn,1)
        conn = fes.conn[i,:];
        econn=zeros(FInt,1,nedges);
        for J = 1:nedges
            ev=conn[ec[J,:]];
            h,n=MeshUtilModule.findhyperface!(edges, ev);
            econn[J]=n;
        end
        inn=size(xyz,1)-size(fes.conn,1)+i
        xyz[inn,:]=mean(xyz[conn[:],:],1); # interior node
        #h,inn=MeshUtilModule.findhyperface!(faces, conn);
        nconn[nc,:] =[conn[1] econn[1] inn econn[4]];
        nc= nc+ 1;
        nconn[nc,:] =[conn[2] econn[2] inn econn[1]];
        nc= nc+ 1;
        nconn[nc,:] =[conn[3] econn[3] inn econn[2]];
        nc= nc+ 1;
        nconn[nc,:] =[conn[4] econn[4] inn econn[3]];
        nc= nc+ 1;
    end
    fens =FENodeSetModule.FENodeSet(xyz);
    nfes = FESetModule.FESetQ4(nconn);
    return fens,nfes            # I think I should not be overwriting the input!
end

"""
    Q8annulus(rin::FFlt, rex::FFlt, nr::FInt, nc::FInt, Angl::FFlt)

Mesh of an annulus segment.

Mesh of an annulus segment, centered at the origin, with internal radius
rin, and  external radius rex, and  development angle Angl. Divided into
elements: nr, nc in the radial and circumferential direction
respectively.
"""
function Q8annulus(rin::FFlt, rex::FFlt, nr::FInt, nc::FInt, Angl::FFlt)
    trin=min(rin,rex);
    trex=max(rin,rex);
    fens,fes = Q8block(trex-trin,Angl,nr,nc);
    xy=fens.xyz;
    for i=1:FENodeSetModule.count(fens)
        r=trin+xy[i,1]; a=xy[i,2];
        xy[i,:]=[r*cos(a) r*sin(a)];
    end
    fens.xyz=xy;
    return fens,fes
end

end
