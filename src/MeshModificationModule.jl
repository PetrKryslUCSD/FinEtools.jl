"""
    MeshModificationModule

Module for mesh modification operations.
"""
module MeshModificationModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FESetModule: FESet, count, boundaryconn, boundaryfe, updateconn!, connasarray, fromarray!
import FinEtools.FENodeSetModule: FENodeSet
import FinEtools.BoxModule: boundingbox, inflatebox!, intersectboxes, inbox
import FinEtools.MeshSelectionModule: connectednodes
using Base.Sort
using Base.Order
import LinearAlgebra: norm, svd, dot, eigen
import Random: randperm

"""
    interior2boundary(interiorconn::Array{Int, 2}, extractb::Array{Int, 2})

Extract the boundary connectivity from the connectivity of the interior.
"""
function interior2boundary(interiorconn::Array{Int, 2}, extractb::Array{Int, 2})
    hypf = interiorconn[:, extractb[1, :]]
    for i = 2:size(extractb, 1)
        hypf = vcat(hypf, interiorconn[:, extractb[i, :]])
    end
    return _myunique2(hypf);
end

"""
meshboundary(fes::T) where {T<:FESet}

Extract the boundary finite elements from a mesh.

Extract the finite elements of manifold dimension (n-1) from the
supplied finite element set of manifold dimension (n).
"""
function meshboundary(fes::T) where {T<:FESet}
    # Form all hyperfaces, non-duplicates are boundary cells
    hypf = boundaryconn(fes);    # get the connectivity of the boundary elements
    bdryconn = _myunique2(hypf);
    make = boundaryfe(fes);     # get the function that can make a boundary element
    return make(bdryconn);
end

function _mysortrows(A::FIntMat)
    # Sort the rows of A by sorting each column from back to front.

    m,n = size(A);

    indx =  zeros(FInt,m); sindx = zeros(FInt,m)
    for i=1:m
        indx[i]=i
    end
    nindx =  zeros(FInt,m);
    col = zeros(FInt,m)
    for c = n:-1:1
        for i=1:m
            col[i]=A[indx[i],c]
        end
        #Sorting a column vector is much faster than sorting a column matrix
        sindx=sortperm(col,alg=QuickSort);
        #sortperm!(sindx,col,alg=QuickSort); # available for 0.4, slightly faster
        #indx=indx[sindx] # saving allocations by using the below loops
        for i=1:m
            nindx[i]=indx[sindx[i]]
        end
        for i=1:m
            indx[i]=nindx[i]
        end
    end

    return A[indx,:]
end

function _mysortdim2!(A::FIntMat)
    # Sort each row  of A in ascending order.

    m,n = size(A);
    r = zeros(FInt,n)
   @inbounds for k = 1:m
        for i=1:n
            r[i]=A[k,i]
        end
        sort!(r);
        for i=1:n
            A[k,i]=r[i]
        end
    end
    return A
end

function  _myunique2(A::FIntVec)
    return _myunique2(reshape(A, length(A), 1))
end

function  _myunique2(A::FIntMat) # speeded up; now the bottleneck is _mysortrows
    #println("size(A)=$(size(A))")
    maxA=maximum(A[:])::FInt
    sA=deepcopy(A)
    #@time
    sA=_mysortdim2!(sA)::FIntMat;#this is fast
    #@time sA=sort(A,2,alg=QuickSort)::FIntMat;#this is slow
    sA= [sA broadcast(+, 1:size(A,1), maxA)]::FIntMat
    #@time
    sA =_mysortrows(sA); # this now takes the majority of time, but much less than the function below
    #@time sA  = sortrows(sA,alg=QuickSort);;#this is slow
    rix=sA[:,end];
    broadcast!(-, rix, rix, maxA)
    sA=sA[:,1:end-1];
    d=falses(size(sA,1)-1)
    for k=1:length(d)
        for m=1:size(sA,2)
            if sA[k,m]!=sA[k+1,m]
                d[k]=true;
                break;
            end
        end
    end
    #d=(sA[1:end-1,:].!=sA[2:end,:]); # element-wise comparison!
    ad=zeros(FInt,size(d,1)+1)
    ad[1]=1;
    for k=2:length(ad)
        for m=1:size(d,2)
            if d[k-1,m]!=0
                ad[k]=1;
                break;
            end
        end
    end
    #ad=map((x) -> (x?1:0),[true; any(d,2)]);
    iu=trues(length(ad))
    for k=1:(length(ad)-1)
        ad[k]=ad[k]+ad[k+1]
        iu[k]=(ad[k]>1)
    end
    ad[end]=ad[end]+1;
    iu[end]=(ad[end]>1)
    #iu =map((x) -> (x>1? true: false),(ad + [ad[2:end];1]));
    Out =A[rix[iu],:];
    return Out
end

# ### This code is correct, but very slow.
# function  _myunique1(A::FIntMat)
#     maxA=maximum(A[:])
#     sA=sort(A,2);# most time spent here
#     sA= [sA (1:size(A,1))+maxA]
#     sA  = sortrows(sA);;#and here
#     rix=sA[:,end]; rix=rix[:]-maxA;
#     sA=sA[:,1:end-1];
#     d=(sA[1:end-1,:].!=sA[2:end,:]); # element-wise comparison!
#     ad=map((x) -> (x ? 1 : 0),[true; any(d,2)]);
#     iu =map((x) -> (x>1 ? true: false),(ad + [ad[2:end];1]));
#     Out =A[rix[iu[:]],:];
#     return Out
# end

"""
    fusenodes(fens1::FENodeSet, fens2::FENodeSet, tolerance:: FFlt)

Fuse together nodes from two node sets.

Fuse two node sets. If necessary, by gluing together nodes located within
tolerance of each other. The two node sets, `fens1` and `fens2`,  are fused
together by merging the nodes that fall within a box of size `tolerance`. The
merged node set, `fens`, and the new  indexes of the nodes in the set `fens1`
are returned.

The set `fens2` will be included unchanged, in the same order,
in the node set `fens`.
The indexes of the node set `fens1` will have changed.

### Example:
After the call to this function we have
`k=new_indexes_of_fens1_nodes[j]` is the node in the node set fens which
used to be node `j` in node set `fens1`.
The finite element set connectivity that used to refer to `fens1`
needs to be updated to refer to the same nodes in  the set `fens` as
     `updateconn!(fes, new_indexes_of_fens1_nodes);`
"""
function fusenodes(fens1::FENodeSet, fens2::FENodeSet, tolerance:: FFlt)
    @assert size(fens1.xyz, 2) == size(fens2.xyz, 2)
    dim::FInt = size(fens1.xyz,2);
    nn1::FInt = count(fens1)
    nn2::FInt = count(fens2)
    xyz1 = zeros(FFlt,nn1,dim); copyto!(xyz1, fens1.xyz)#::FFltMat = copy(fens1.xyz::FFltMat)
    id1 = collect(1:nn1);
    xyz2 = zeros(FFlt,nn2,dim); copyto!(xyz2, fens2.xyz)#xyz2::FFltMat = copy(fens2.xyz::FFltMat)
    id2 = collect(1:nn2);
    # Decide which nodes should be checked for proximity
    ib::FFltVec = intersectboxes(inflatebox!(boundingbox(xyz1), tolerance), inflatebox!(boundingbox(xyz2), tolerance))
    node1in = fill(false, nn1);
    node2in = fill(false, nn2);
    if length(ib) > 0
        for i=1:nn1
            node1in[i] = inbox(ib, @view xyz1[i, :])
        end
        for i=1:nn2
            node2in[i] = inbox(ib, @view xyz2[i, :])
        end
    end 
    # Mark nodes from the first array that are duplicated in the second
    if (tolerance > 0.0) # should we attempt to merge nodes?
        for i=1:nn1
            if node1in[i]
                breakoff = false
                for rx=1:nn2
                    if node2in[rx]
                        distance::FFlt= 0.0
                        for cx=1:dim
                            distance = distance + abs(xyz2[rx,cx]-xyz1[i,cx]);
                            if (distance >= tolerance) # shortcut: if the distance is already too large, stop checking
                                break
                            end
                        end
                        if (distance < tolerance)
                            id1[i] = -rx; breakoff = true;
                        end
                    end 
                    if breakoff
                        break
                    end
                end
            end
        end
    end
    # Generate  fused arrays of the nodes. First copy in the nodes from the second set...
    xyzm = zeros(FFlt,nn1+nn2,dim);
    for rx = 1:nn2
        for cx = 1:dim
            xyzm[rx,cx] = xyz2[rx,cx];
        end
    end
    idm = zeros(FInt,nn1+nn2);
    for rx = 1:nn2
        idm[rx] = rx;
    end
    mid=nn2+1;
    # ...and then we add in only non-duplicated nodes from the first set
    for i=1:nn1 
        if id1[i]>0
            id1[i] = mid;
            idm[mid] = mid;
            for cx = 1:dim
                xyzm[mid,cx] = xyz1[i,cx];
            end
            mid = mid+1;
        else
            id1[i] = id2[-id1[i]];
        end
    end
    nnodes = mid-1;
    xyzm = xyzm[1:nnodes,:];

    # Create the fused Node set
    fens = FENodeSet(xyzm);
    # The Node set 1 numbering will change
    new_indexes_of_fens1_nodes = id1[:];
    # The node set 2 numbering stays the same
    return fens, new_indexes_of_fens1_nodes
end

"""
    compactnodes(fens::FENodeSet, connected::Vector{Bool})

Compact the finite element node set by deleting unconnected nodes.

`fens` = array of finite element nodes
`connected` = The array element `connected[j]` is either 0 (when `j` is an
  unconnected node), or a positive number (when node `j` is connected to
  other nodes by at least one finite element)

### Output:
`fens` = new set of finite element nodes
`new_numbering`= array which tells where in the new `fens` array the
     connected nodes are (or 0 when the node was unconnected). For instance,
     node 5 was connected, and in the new array it is the third node: then
     `new_numbering[5]` is 3.

### Examples:
Let us say there are nodes not connected to any finite element that you
would like to remove from the mesh: here is how that would be
accomplished.
```
connected = findunconnnodes(fens, fes);
fens, new_numbering = compactnodes(fens, connected);
fes = renumberconn!(fes, new_numbering);
```
Finally, check that the mesh is valid:
```
validate_mesh(fens, fes);
```
"""
function compactnodes(fens::FENodeSet, connected::BitArray{1})
    @assert length(connected) == count(fens)
    new_numbering = zeros(FInt,count(fens),1);
    nxyz = deepcopy(fens.xyz);
    id=1;
    for i=1:length(connected)
        if (connected[i])
            new_numbering[i] = id;
            nxyz[id,:] = fens.xyz[i,:];
            id=id+1;
        end
    end
    #new_numbering = new_numbering[1:id-1];
    fens = FENodeSet(nxyz[1:id-1,:]);
    return fens, vec(new_numbering)
end

"""
    mergemeshes(fens1::FENodeSet, fes1::T1,
      fens2::FENodeSet, fes2::T2, tolerance::FFlt) where {T1<:FESet,T2<:FESet}

Merge together two meshes.

Merge two meshes together by gluing together nodes within tolerance. The
two meshes, `fens1`, `fes1`, and `fens2`, `fes2`, are glued together by merging
the nodes that fall within a box of size `tolerance`. If `tolerance` is set
to zero, no merging of nodes is performed; the two meshes are simply
concatenated together.

The merged node set, `fens`, and the two finite element sets with
renumbered  connectivities are returned.

Important notes: On entry into this function the connectivity of `fes1`
point into `fens1` and the connectivity of `fes2` point into `fens2`. After
this function returns the connectivity of both `fes1` and `fes2` point into
`fens`. The order of the nodes of the node set `fens1` in the resulting set
`fens` will have changed, whereas the order of the nodes of the node set
`fens2` is are guaranteed to be the same. Therefore, the connectivity of
`fes2` will in fact remain the same.
"""
function mergemeshes(fens1::FENodeSet, fes1::T1,
    fens2::FENodeSet, fes2::T2, tolerance::FFlt) where {T1<:FESet,T2<:FESet}
    # Fuse the nodes
    # @code_warntype fusenodes(fens1, fens2, tolerance);
    fens, new_indexes_of_fens1_nodes = fusenodes(fens1, fens2, tolerance);
    # Renumber the finite elements
    newfes1 = deepcopy(fes1)
    updateconn!(newfes1, new_indexes_of_fens1_nodes);
    # Note that now the connectivity of both fes1 and fes2 point into
    # fens.
    return fens, newfes1, fes2
end

"""
    mergenmeshes(meshes::Array{Tuple{FENodeSet, FESet}}, tolerance::FFlt)

Merge several meshes together.

The meshes are glued together by merging the nodes that fall within 
a box of size `tolerance`. If `tolerance` is set to zero, no merging of 
nodes is performed; the nodes from the meshes are simply concatenated together.

## Output
The merged node set, `fens`, and an array of finite element sets with
renumbered  connectivities are returned.
"""
function mergenmeshes(meshes::Array{Tuple{FENodeSet, FESet}}, tolerance::FFlt)
    outputfes = Array{FESet,1}()
    if (length(meshes)) == 1 # A single mesh, package output and return
        fens, fes = meshes[1];
        push!(outputfes, fes)
        return fens, outputfes
    end
    # Multiple meshes: process
    fens, fes = meshes[1];
    push!(outputfes, fes)
    for j=2:length(meshes)
        fens1, fes1 = meshes[j];
        fens, new_indexes_of_fens1_nodes = fusenodes(fens1, fens, tolerance);
        updateconn!(fes1,new_indexes_of_fens1_nodes);
        push!(outputfes, fes1)
    end
    return fens, outputfes
end

"""
    mergenodes(fens::FENodeSet, fes::FESet, tolerance::FFlt)

Merge together  nodes of a single node set.

Merge by gluing together nodes from a single node set located within
tolerance of each other. The nodes are glued together by merging the
nodes that fall within a box of size `tolerance`. The merged node
set, `fens`, and the finite element set, `fes`, with renumbered  connectivities
are returned.

Warning: This tends to be an expensive operation!
"""
function mergenodes(fens::FENodeSet, fes::FESet, tolerance::FFlt)
    maxnn = count(fens) + 1
    xyz1 = fens.xyz;
    dim  = size(xyz1,2);
    id1 = collect(1:count(fens));
    d = zeros(size(xyz1,1));
    # Mark nodes from the array that are duplicated
    for i = 1:count(fens)
        if (id1[i] > 0) # This node has not yet been marked for merging
            XYZ = reshape(xyz1[i,:], 1, dim);
            copyto!(d, sum(abs.(xyz1 .- XYZ), dims = 2)); #find the distances along  coordinate directions
            minn = maxnn
            @inbounds for jx = 1:length(d)
                if d[jx] < tolerance
                    minn = min(jx, minn);
                    id1[jx] = -minn;
                    id1[minn] = minn;
                end
            end
        end
    end
    # Generate  merged arrays of the nodes
    xyzm = zeros(FFlt,count(fens),dim);
    mid = 1;
    for i = 1:count(fens) # and then we pick only non-duplicated fens1
        if id1[i] > 0 # this node is the master
            id1[i] = mid;
            xyzm[mid,:] = xyz1[i,:];
            mid = mid+1;
        else # this node is the slave
            id1[i] = id1[-id1[i]];
        end
    end
    nnodes = mid-1;
    xyzm = xyzm[1:nnodes,:];
    # Renumber the cells
    conns = connasarray(fes);
    for i = 1:size(conns,1)
        conns[i,:] = id1[conns[i,:]];
    end
    fes = fromarray!(fes, conns)

    fens = FENodeSet(xyzm[1:nnodes,:]);

    return fens,fes
end

"""
    mergenodes(fens::FENodeSet, fes::FESet, tolerance::FFlt, candidates::FIntVec)

Merge together  nodes of a single node set.

Similar to `mergenodes(fens::FENodeSet, fes::FESet, tolerance::FFlt)`, but only
the candidate nodes are considered for merging. This can potentially speed up the operation by orders of magnitude.
"""
function mergenodes(fens::FENodeSet, fes::FESet, tolerance::FFlt, candidates::FIntVec)
    maxnn = count(fens) + 1
    xyz1 = fens.xyz;
    dim  = size(xyz1,2);
    id1 = collect(1:count(fens));
    d = fill(100.0*tolerance, size(xyz1, 1))
    # Mark nodes from the array that are duplicated
    for ic = 1:length(candidates)
        i = candidates[ic]
        if (id1[i] > 0) # This node has not yet been marked for merging
            XYZ = xyz1[i,:];
            minn = maxnn
            for kx = candidates
                d[kx] = sum(abs.(xyz1[kx,:] .- XYZ))
            end
            @inbounds for jx = candidates
                if d[jx] < tolerance
                    minn = min(jx, minn);
                    id1[jx] = -minn;
                    id1[minn] = minn;
                end
            end
        end
    end
    # Generate  merged arrays of the nodes
    xyzm = zeros(FFlt,count(fens),dim);
    mid = 1;
    for i = 1:count(fens) # and then we pick only non-duplicated fens1
        if id1[i] > 0 # this node is the master
            id1[i] = mid;
            xyzm[mid,:] = xyz1[i,:];
            mid = mid+1;
        else # this node is the slave
            id1[i] = id1[-id1[i]];
        end
    end
    nnodes = mid-1;
    xyzm = xyzm[1:nnodes,:];
    # Renumber the cells
    conns = connasarray(fes);
    for i = 1:count(fes)
        conn = conns[i,:];
        conns[i,:] = id1[conn];
    end
    fes = fromarray!(fes, conns);

    fens = FENodeSet(xyzm[1:nnodes,:]);

    return fens,fes
end

"""
    renumberconn!(fes::FESet, new_numbering::FIntVec)

Renumber the nodes in the connectivity of the finite elements based on a new
numbering for the nodes.

`fes` =finite element set
`new_numbering` = new serial numbers for the nodes.  The connectivity
          should be changed as `conn[j]` --> `new_numbering[conn[j]]`

Let us say there are nodes not connected to any finite element that you would
like to remove from the mesh: here is how that would be accomplished.
```
connected = findunconnnodes(fens, fes);
fens, new_numbering = compactfens(fens, connected);
fes = renumberconn!(fes, new_numbering);
```
Finally, check that the mesh is valid:
```julia
validate_mesh(fens, fes);
```
"""
function renumberconn!(fes::FESet, new_numbering::FIntVec)
    conn = connasarray(fes)
    for i=1:size(conn,1)
        c = conn[i,:];
        conn[i,:] = new_numbering[c];
    end
    return fromarray!(fes, conn)
end

"""
    vsmoothing(v::FFltMat, t::FIntMat; options...)

Internal routine for mesh smoothing.

Keyword options:
`method` = :taubin (default) or :laplace
`fixedv` = Boolean array, one entry per vertex: is the vertex immovable (true)
    or movable  (false)
`npass` = number of passes (default 2)
"""
function vsmoothing(v::FFltMat, t::FIntMat; kwargs...)
    fixedv = falses(size(v,1))
    npass = 2;
    method =:taubin;
    for apair in pairs(kwargs)
        sy, val = apair
        if sy==:method
            method = val
        elseif sy==:fixedv
            fixedv .= val
        elseif sy==:npass
            npass = val
        end
    end

    nv = deepcopy(v)
    # find neighbors for the given connections
    vneigh =  vertexneighbors(t,size(v,1));
    # Smoothing considering all connections through the volume
    if (method == :taubin)
        nv =  smoothertaubin(v,vneigh,fixedv,npass,0.5,-0.5);
    elseif (method == :laplace)
        nv =  smootherlaplace(v,vneigh,fixedv,npass,0.5,-0.5);
    end
    # return new vertex locations
    return nv
end

"""
    meshsmoothing(fens::FENodeSet, fes::T; options...) where {T<:FESet}

General smoothing of meshes.

## Keyword options:
`method` = :taubin (default) or :laplace
`fixedv` = Boolean array, one entry per vertex: is the vertex immovable (true)
    or movable  (false)
`npass` = number of passes (default 2)

## Return
The modified  node set.
"""
function meshsmoothing(fens::FENodeSet, fes::T; options...) where {T<:FESet}
    v = deepcopy(fens.xyz)
    v = vsmoothing(v, connasarray(fes); options...)
    copyto!(fens.xyz, v)
    return fens
end

function  smoothertaubin(vinp::FFltMat, vneigh::Array{FIntVec,1}, fixedv::T, npass::FInt, lambda::FFlt, mu::FFlt) where {T}
    v=deepcopy(vinp);
    nv=deepcopy(vinp);
    for I= 1:npass
        o=randperm(length(vneigh));
        damping_factor=lambda;
        for k= 1:length(vneigh)
            r=o[k];
            n=vneigh[r];
            if (length(n)>1) && (!fixedv[r])
                ln1 = (length(n)-1)
                nv[r,:] .= (1-damping_factor)*vec(v[r,:]) + damping_factor*(vec(sum(v[n,:], dims = 1)) - vec(v[r,:]))/ln1;
            end
        end
        v=deepcopy(nv);
        damping_factor=mu;
        for k= 1:length(vneigh)
            r=o[k];
            n=vneigh[r];
            if (length(n)>1) && (!fixedv[r])
                ln1 = (length(n)-1)
                nv[r,:] .= (1-damping_factor)*vec(v[r,:]) + damping_factor*(vec(sum(v[n,:], dims = 1)) - vec(v[r,:]))/ln1;
            end
        end
        v=deepcopy(nv);
    end
    return nv
end

function   smootherlaplace(vinp::FFltMat, vneigh::Array{FIntVec,1}, fixedv::T, npass::FInt, lambda::FFlt,mu::FFlt) where {T}
    v=deepcopy(vinp);
    nv=deepcopy(vinp);
    damping_factor=lambda;
    for I= 1:npass
        o=randperm(length(vneigh));
        for k= 1:length(vneigh)
            r=o[k];
            n=vneigh[r];
            if (length(n)>1) && (!fixedv[r])
                ln1 = (length(n)-1)
                nv[r,:] = (1-damping_factor)*vec(v[r,:]) + damping_factor*(vec(sum(v[n,:], dims = 1))-vec(v[r,:]))/ln1;
            end
        end
        v=deepcopy(nv);
    end
    return nv
end

"""
    vertexneighbors(conn::FIntMat, nvertices::FInt)

Find the node neighbors in the mesh. 

Return an array of integer vectors, element I holds an array of numbers of nodes
which are connected to node I (including node I).  
"""
function vertexneighbors(conn::FIntMat, nvertices::FInt)
    vn = FIntVec[]; sizehint!(vn, nvertices)
    for I= 1:nvertices
        push!(vn, FInt[]);          # preallocate
    end
    for I= 1:size(conn,1)
        for r= 1:size(conn,2)
            append!(vn[conn[I,r]],vec(conn[I,:]));
        end
    end
    for I= 1:length(vn)
        vn[I]=unique(vn[I]);
    end
    return vn
end

"""
    mirrormesh(fens::FENodeSet, fes::T, Normal::FFltVec,
      Point::FFltVec; kwargs...) where {T<:FESet}

Mirror a 2-D mesh in a plane given by its normal and one point.

Warning: The code to relies on the numbering of the cells: to reverse the
orientation of the mirrored cells, the connectivity is listed in reverse
order.   If the mirrored cells do not follow this rule (for instance hexahedra
for quadrilaterals), their areas/volumes will come out negative. In such a
case the renumbering function of the connectivity needs to be supplied.

For instance: H8 elements require  the renumbering function to be supplied as
fens1,gcells1 = mirror_mesh(fens, gcells,...
          [-1,0,0], [0,0,0], @(c)c([1, 4, 3, 2, 5, 8, 7, 6]));
"""
function mirrormesh(fens::FENodeSet, fes::T, Normal::FFltVec,
    Point::FFltVec; kwargs...) where {T<:FESet}
    # Treat optional arguments.
    # Simply switch the order of nodes.  Works for simplexes...
    renumb(conn) = conn[end:-1:1];
    for apair in pairs(kwargs)
        sy, val = apair
        if sy == :renumb
            renumb = val
        end
    end
    # Make sure we're using a unit normal
    Normal = Normal/norm(Normal);
    Normal = vec(Normal)
    # The point needs to be a row  matrix
    Point = vec(Point)

    fens1 = deepcopy(fens); # the mirrored mesh nodes
    for i = 1:count(fens1)
        a = fens1.xyz[i,:]
        d = dot(vec(a-Point), Normal);
        fens1.xyz[i,:] = a-2*d*Normal;
    end
    # Reconnect the cells
    fes1=deepcopy(fes);
    conn = connasarray(fes1)
    for i=1:size(conn, 1)
        conn[i,:]=renumb(conn[i,:]);
    end
    return fens1, fromarray!(fes1, conn)
end
 
function _nodepartitioning3(fens::FENodeSet, nincluded::Vector{Bool}, npartitions::Int = 2)
    function inertialcutpartitioning!(partitioning, parts, X)
        nspdim = 3
        StaticMoments = fill(zero(FFlt), nspdim, length(parts));
        npart = fill(0, length(parts))
        for spdim = 1:nspdim
            @inbounds for j = 1:size(X, 1)
                if nincluded[j] # Is the node to be included in the partitioning?
                    jp = partitioning[j]
                    StaticMoments[spdim, jp] += X[j, spdim]
                    npart[jp] += 1 # count the nodes in the current partition
                end 
            end
        end
        CG = fill(zero(FFlt), nspdim, length(parts));
        for p = parts
            npart[p] = Int(npart[p] / nspdim)
            CG[:, p] = StaticMoments[:, p] / npart[p] # center of gravity of each partition
        end 
        MatrixMomentOfInertia = fill(zero(FFlt), nspdim, nspdim, length(parts))
        @inbounds for j = 1:size(X, 1)
            if nincluded[j] # Is the node to be included in the partitioning?
                jp = partitioning[j]
                xj, yj, zj = X[j, 1] - CG[1, jp], X[j, 2] - CG[2, jp], X[j, 3] - CG[3, jp]
                MatrixMomentOfInertia[1, 1, jp] += yj^2 + zj^2
                MatrixMomentOfInertia[2, 2, jp] += xj^2 + zj^2
                MatrixMomentOfInertia[3, 3, jp] += yj^2 + xj^2
                MatrixMomentOfInertia[1, 2, jp] -= xj * yj
                MatrixMomentOfInertia[1, 3, jp] -= xj * zj
                MatrixMomentOfInertia[2, 3, jp] -= yj * zj
            end 
        end
        for p = parts
            MatrixMomentOfInertia[2, 1, p] = MatrixMomentOfInertia[1, 2, p]
            MatrixMomentOfInertia[3, 1, p] = MatrixMomentOfInertia[3, 1, p]
            MatrixMomentOfInertia[3, 2, p] = MatrixMomentOfInertia[3, 2, p]
        end 
        longdir = fill(zero(FFlt), nspdim, length(parts))
        for p = parts
            F = eigen(MatrixMomentOfInertia[:, :, p])
            six = sortperm(F.values)
            longdir[:, p] = F.vectors[:, six[1]]
        end 
        toggle = fill(one(FFlt), length(parts)); 
        @inbounds for j = 1:size(X, 1)
            if nincluded[j] # Is the node to be included in the partitioning?
                jp = partitioning[j]
                vx, vy, vz = longdir[:, jp]
                xj, yj, zj = X[j, 1] - CG[1, jp], X[j, 2] - CG[2, jp], X[j, 3] - CG[3, jp]
                d = xj * vx + yj * vy + zj * vz
                c = 0
                if d < 0.0
                    c = 1
                elseif d > 0.0
                    c = 0
                else # disambiguate d[ixxxx] == 0.0
                    c = (toggle[jp] > 0) ? 1 : 0
                    toggle[jp] = -toggle[jp]
                end
                partitioning[j] = 2 * jp - c
            end 
        end 
    end
    
    nlevels = Int(round(ceil(log(npartitions)/log(2))))
    partitioning = fill(1, count(fens))  # start with nodes assigned to partition 1
    for level = 0:1:(nlevels - 1)
        inertialcutpartitioning!(partitioning, collect(1:2^level), fens.xyz)
    end
    return partitioning
end

function _nodepartitioning2(fens::FENodeSet, nincluded::Vector{Bool}, npartitions::Int = 2)
    function inertialcutpartitioning!(partitions, parts, X)
        nspdim = 2
        StaticMoments = fill(zero(FFlt), nspdim, length(parts));
        npart = fill(0, length(parts))
        for spdim = 1:nspdim
            @inbounds for j = 1:size(X, 1)
                if nincluded[j] # Is the node to be included in the partitioning?
                    jp = partitions[j]
                    StaticMoments[spdim, jp] += X[j, spdim]
                    npart[jp] += 1 # count the nodes in the current partition
                end 
            end
        end
        CG = fill(zero(FFlt), nspdim, length(parts));
        for p = parts
            npart[p] = Int(npart[p] / nspdim)
            CG[:, p] = StaticMoments[:, p] / npart[p] # center of gravity of each partition
        end 
        MatrixMomentOfInertia = fill(zero(FFlt), nspdim, nspdim, length(parts))
        @inbounds for j = 1:size(X, 1)
            if nincluded[j] # Is the node to be included in the partitioning?
                jp = partitions[j]
                xj, yj = X[j, 1] - CG[1, jp], X[j, 2] - CG[2, jp]
                MatrixMomentOfInertia[1, 1, jp] += yj^2
                MatrixMomentOfInertia[2, 2, jp] += xj^2
                MatrixMomentOfInertia[1, 2, jp] -= xj * yj
            end 
        end
        for p = parts
            MatrixMomentOfInertia[2, 1, p] = MatrixMomentOfInertia[1, 2, p]
        end 
        longdir = fill(zero(FFlt), nspdim, length(parts))
        for p = parts
            F = eigen(MatrixMomentOfInertia[:, :, p])
            six = sortperm(F.values)
            longdir[:, p] = F.vectors[:, six[1]]
        end 
        toggle = fill(one(FFlt), length(parts)); 
        @inbounds for j = 1:size(X, 1)
            if nincluded[j] # Is the node to be included in the partitioning?
                jp = partitions[j]
                vx, vy = longdir[:, jp]
                xj, yj = X[j, 1] - CG[1, jp], X[j, 2] - CG[2, jp]
                d = xj * vx + yj * vy
                c = 0
                if d < 0.0
                    c = 1
                elseif d > 0.0
                    c = 0
                else # disambiguate d[ixxxx] == 0.0
                    c = (toggle[jp] > 0) ? 1 : 0
                    toggle[jp] = -toggle[jp]
                end
                partitions[j] = 2 * jp - c
            end 
        end 
    end
    
    nlevels = Int(round(ceil(log(npartitions)/log(2))))
    partitions = fill(1, count(fens))  # start with nodes assigned to partition 1
    for level = 0:1:(nlevels - 1)
        inertialcutpartitioning!(partitions, collect(1:2^level), fens.xyz)
    end
    return partitions
end

"""
    nodepartitioning(fens::FENodeSet, nincluded::Vector{Bool}, npartitions)

Compute the inertial-cut partitioning of the nodes.

`nincluded` = Boolean array: is the node to be included in the partitioning or
    not?
`npartitions` = number of partitions, but note that the actual number of
    partitions is going to be a power of two.

The partitioning can be visualized for instance as:
```julia
partitioning = nodepartitioning(fens, npartitions)
partitionnumbers = unique(partitioning)
for gp = partitionnumbers
  groupnodes = findall(k -> k == gp, partitioning)
  File =  "partition-nodes-Dollar(gp).vtk"
  vtkexportmesh(File, fens, FESetP1(reshape(groupnodes, length(groupnodes), 1)))
end 
File =  "partition-mesh.vtk"
vtkexportmesh(File, fens, fes)
@async run(`"paraview.exe" DollarFile`)
```
"""
function nodepartitioning(fens::FENodeSet, nincluded::Vector{Bool}, npartitions::Int = 2)
    @assert npartitions >= 2
    if size(fens.xyz, 2) == 3
        return _nodepartitioning3(fens, nincluded, npartitions)
    elseif size(fens.xyz, 2) == 2
        return _nodepartitioning2(fens, nincluded, npartitions)
    else 
        @warn "Not implemented for 1D"
    end
end 

"""
    nodepartitioning(fens::FENodeSet, npartitions)

Compute the inertial-cut partitioning of the nodes.

`npartitions` = number of partitions, but note that the actual number of
partitions is going to be a power of two.

In this variant all the nodes are to be included in the partitioning.

The partitioning can be visualized for instance as:
```julia
partitioning = nodepartitioning(fens, npartitions)
partitionnumbers = unique(partitioning)
for gp = partitionnumbers
  groupnodes = findall(k -> k == gp, partitioning)
  File =  "partition-nodes-Dollar(gp).vtk"
  vtkexportmesh(File, fens, FESetP1(reshape(groupnodes, length(groupnodes), 1)))
end 
File =  "partition-mesh.vtk"
vtkexportmesh(File, fens, fes)
@async run(`"paraview.exe" DollarFile`)
```
"""
function nodepartitioning(fens::FENodeSet, npartitions::Int = 2)
    @assert npartitions >= 2
    nincluded = fill(true, count(fens)) # The default is all nodes are included in the partitioning.
    if size(fens.xyz, 2) == 3
        return _nodepartitioning3(fens, nincluded, npartitions)
    elseif size(fens.xyz, 2) == 2
        return _nodepartitioning2(fens, nincluded, npartitions)
    else 
        @warn "Not implemented for 1D"
    end
end 

"""
    nodepartitioning(fens::FENodeSet, fesarr, npartitions::Vector{Int})

Compute the inertial-cut partitioning of the nodes.

`fesarr` = array of finite element sets that represent regions
`npartitions` = array of the number of partitions in each region. However 
note that the actual number of partitions is going to be a power of two.

The partitioning itself is carried out by `nodepartitioning()` with 
a list of nodes to be included in the partitioning. For each region I
the nodes included in the partitioning are those connected to 
the elements of that region, but not to elements that belong to 
any of the previous regions, 1…I-1.
"""
function nodepartitioning(fens::FENodeSet, fesarr, npartitions::Vector{Int})
    @assert length(fesarr) == length(npartitions)
    # Partitioning of all the nodes
    partitioning = fill(0, count(fens))
    # Find the partitioning of the nodes in FESet 1
    nincludedp = fill(false, count(fens))
    for i = connectednodes(fesarr[1]) # For nodes connected by region 1
        nincludedp[i] = true
    end
    partitioning1 = nodepartitioning(fens, nincludedp, npartitions[1])
    totnpartitions = maximum(partitioning1)
    # Transfer the partitioning of region 1 into the overall partitioning
    partitioning[nincludedp] = partitioning1[nincludedp]
    for i = 2:length(npartitions)
        # Find the partitioning of the nodes in FESet i, but not in the preceding sets
        nincluded = fill(false, count(fens))
        for j = connectednodes(fesarr[i]) # For nodes connected by region i
            nincluded[j] = !nincludedp[j] # Not included previously
        end
        partitioning1 = nodepartitioning(fens, nincluded, npartitions[i])
        np = maximum(partitioning1) # Number of partitions for region i
        partitioning1 = partitioning1 .+ totnpartitions # shift by the number of partitions in previous regions
        totnpartitions = totnpartitions + np # Increment the number of partitions
        partitioning[nincluded] = partitioning1[nincluded] # Update overall partitioning 
        @. nincludedp = nincluded | nincludedp  # Update the list of nodes that have already been included
    end
    return partitioning
end 

end
