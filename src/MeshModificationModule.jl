"""
    MeshModificationModule

Module for mesh modification operations.
"""
module MeshModificationModule

__precompile__(true)

import ..FESetModule:
    AbstractFESet,
    count,
    boundaryconn,
    boundaryfe,
    updateconn!,
    connasarray,
    fromarray!,
    centroidparametric,
    bfun,
    bfundpar,
    subset
import ..FENodeSetModule: FENodeSet
import ..BoxModule: boundingbox, inflatebox!, intersectboxes, inbox
import ..MeshSelectionModule: connectednodes, selectelem, findunconnnodes
import ..SurfaceNormalModule: SurfaceNormal, updatenormal!
using Base.Sort
using Base.Order
import LinearAlgebra: norm, svd, dot, eigen
import Random: randperm
using Statistics: mean
using SparseArrays

"""
    interior2boundary(interiorconn::Array{IT,2}, extractb::Array{IT,2}) where {IT<:Integer}

Extract the boundary connectivity from the connectivity of the interior.

`extractb` = array that defines in which order the bounding faces are traversed. For
    example, for tetrahedra this array is
    ```
    extractb = [1 3 2; 1 2 4; 2 3 4; 1 4 3]
    ```
"""
function interior2boundary(
    interiorconn::Array{IT,2},
    extractb::Array{IT,2},
) where {IT<:Integer}
    hypf = interiorconn[:, extractb[1, :]]
    for i = 2:size(extractb, 1)
        hypf = vcat(hypf, interiorconn[:, extractb[i, :]])
    end
    return _myunique2(hypf)
end

"""
    meshboundary(fes::T) where {T<:AbstractFESet}

Compute the boundary of a mesh defined by the given finite element set.

# Arguments
- `fes::T`: The finite element set representing the mesh.

# Returns
The boundary of the mesh.

Extract the finite elements of manifold dimension (n-1) from the
supplied finite element set of manifold dimension (n).
"""
function meshboundary(fes::T) where {T<:AbstractFESet}
    # Form all hyperfaces, non-duplicates are boundary cells
    hypf = boundaryconn(fes)    # get the connectivity of the boundary elements
    bdryconn = _myunique2(hypf)
    make = boundaryfe(fes)     # get the function that can make a boundary element
    return make(bdryconn)
end

function _mysortrows(A::Matrix{IT}) where {IT<:Integer}
    # Sort the rows of A by sorting each column from back to front.

    m, n = size(A)

    indx = zeros(IT, m)
    sindx = zeros(IT, m)
    for i = 1:m
        indx[i] = i
    end
    nindx = zeros(IT, m)
    col = zeros(IT, m)
    for c = n:-1:1
        for i = 1:m
            col[i] = A[indx[i], c]
        end
        #Sorting a column vector is much faster than sorting a column matrix
        sindx = sortperm(col, alg = QuickSort)
        #sortperm!(sindx,col,alg=QuickSort); # available for 0.4, slightly faster
        #indx=indx[sindx] # saving allocations by using the below loops
        for i = 1:m
            nindx[i] = indx[sindx[i]]
        end
        for i = 1:m
            indx[i] = nindx[i]
        end
    end

    return A[indx, :]
end

function _mysortdim2!(A::Matrix{IT}) where {IT<:Integer}
    # Sort each row  of A in ascending order.

    m, n = size(A)

    r = zeros(IT, n)
    @inbounds for k = 1:m
        for i = 1:n
            r[i] = A[k, i]
        end
        sort!(r)
        for i = 1:n
            A[k, i] = r[i]
        end
    end
    return A
end

function _myunique2(A::Vector{IT}) where {IT<:Integer}
    return _myunique2(reshape(A, length(A), 1))
end

function _myunique2(A::Matrix{IT}) where {IT<:Integer} # speeded up; now the bottleneck is _mysortrows
    Out = A[_myunique2index(A), :]
end

function _myunique2index(A::Matrix{IT}) where {IT<:Integer} # speeded up; now the bottleneck is _mysortrows
    #println("size(A)=$(size(A))")
    maxA = maximum(A[:])
    sA = deepcopy(A)
    #@time
    sA = _mysortdim2!(sA) #this is fast
    #@time sA=sort(A,2,alg=QuickSort)::FIntMat;#this is slow
    sA = [sA broadcast(+, 1:size(A, 1), maxA)]
    #@time
    sA = _mysortrows(sA) # this now takes the majority of time, but much less than the function below
    #@time sA  = sortrows(sA,alg=QuickSort);;#this is slow
    rix = sA[:, end]
    broadcast!(-, rix, rix, maxA)
    sA = sA[:, 1:(end-1)]
    d = falses(size(sA, 1) - 1)
    for k in eachindex(d)
        for m in axes(sA, 2)
            if sA[k, m] != sA[k+1, m]
                d[k] = true
                break
            end
        end
    end
    #d=(sA[1:end-1,:].!=sA[2:end,:]); # element-wise comparison!
    ad = zeros(IT, size(d, 1) + 1)
    ad[1] = 1
    for k = 2:lastindex(ad)
        for m in axes(d, 2)
            if d[k-1, m] != 0
                ad[k] = 1
                break
            end
        end
    end
    #ad=map((x) -> (x?1:0),[true; any(d,2)]);
    iu = trues(length(ad))
    for k = 1:(lastindex(ad)-1)
        ad[k] = ad[k] + ad[k+1]
        iu[k] = (ad[k] > 1)
    end
    ad[end] = ad[end] + 1
    iu[end] = (ad[end] > 1)
    #iu =map((x) -> (x>1? true: false),(ad + [ad[2:end];1]));
    return rix[iu]
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
    fusenodes(fens1::FENodeSet{T}, fens2::FENodeSet{T}, tolerance::T) where {T<:Number}

Fuse together nodes from two node sets.

Fuse two node sets. If necessary, by gluing together nodes located within
tolerance of each other. The two node sets, `fens1` and `fens2`,  are fused
together by merging the nodes that fall within a box of size `tolerance`. The
merged node set, `fens`, and the new  indexes of the nodes in the set `fens1`
are returned.

The set `fens2` will be included unchanged, in the same order,
in the node set `fens`.
The indexes of the node set `fens1` will have changed.

# Example
After the call to this function we have
`k=new_indexes_of_fens1_nodes[j]` is the node in the node set `fens` which
used to be node `j` in node set `fens1`.
The finite element set connectivity that used to refer to `fens1`
needs to be updated to refer to the same nodes in  the set `fens` as
     `updateconn!(fes, new_indexes_of_fens1_nodes);`
"""
function fusenodes(fens1::FENodeSet{T}, fens2::FENodeSet{T}, tolerance::T) where {T<:Number}
    (size(fens1.xyz, 2) == size(fens2.xyz, 2)) || error("The two node sets must have the same number of space dimensions")
    dim = size(fens1.xyz, 2)
    nn1 = count(fens1)
    nn2 = count(fens2)
    xyz1 = zeros(T, nn1, dim)
    copyto!(xyz1, fens1.xyz)
    id1 = collect(1:nn1)
    xyz2 = zeros(T, nn2, dim)
    copyto!(xyz2, fens2.xyz)
    id2 = collect(1:nn2)
    # Decide which nodes should be checked for proximity
    ib = intersectboxes(
        inflatebox!(boundingbox(xyz1), tolerance),
        inflatebox!(boundingbox(xyz2), tolerance),
    )
    node1in = fill(false, nn1)
    node2in = fill(false, nn2)
    if length(ib) > 0
        for i = 1:nn1
            node1in[i] = inbox(ib, @view xyz1[i, :])
        end
        for i = 1:nn2
            node2in[i] = inbox(ib, @view xyz2[i, :])
        end
    end
    # Mark nodes from the first array that are duplicated in the second
    if (tolerance > 0.0) # should we attempt to merge nodes?
        for i = 1:nn1
            if node1in[i]
                breakoff = false
                for rx = 1:nn2
                    if node2in[rx]
                        distance = T(0.0)
                        for cx = 1:dim
                            distance = distance + abs(xyz2[rx, cx] - xyz1[i, cx])
                            if (distance >= tolerance) # shortcut: if the distance is already too large, stop checking
                                break
                            end
                        end
                        if (distance < tolerance)
                            id1[i] = -rx
                            breakoff = true
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
    xyzm = zeros(eltype(xyz1), nn1 + nn2, dim)
    for rx = 1:nn2
        for cx = 1:dim
            xyzm[rx, cx] = xyz2[rx, cx]
        end
    end
    idm = zeros(Int, nn1 + nn2)
    for rx = 1:nn2
        idm[rx] = rx
    end
    mid = nn2 + 1
    # ...and then we add in only non-duplicated nodes from the first set
    for i = 1:nn1
        if id1[i] > 0
            id1[i] = mid
            idm[mid] = mid
            for cx = 1:dim
                xyzm[mid, cx] = xyz1[i, cx]
            end
            mid = mid + 1
        else
            id1[i] = id2[-id1[i]]
        end
    end
    nnodes = mid - 1
    xyzm = xyzm[1:nnodes, :]

    # Create the fused Node set
    fens = FENodeSet(xyzm)
    # The Node set 1 numbering will change
    new_indexes_of_fens1_nodes = id1[:]
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

# Output
`fens` = new set of finite element nodes
`new_numbering`= array which tells where in the new `fens` array the
     connected nodes are (or 0 when the node was unconnected). For instance,
     node 5 was connected, and in the new array it is the third node: then
     `new_numbering[5]` is 3.

# Examples
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
    (length(connected) == count(fens)) || error("Length of the connected array must be == number of nodes")
    new_numbering = zeros(Int, count(fens), 1)
    nxyz = deepcopy(fens.xyz)
    id = 1
    for i in eachindex(connected)
        if (connected[i])
            new_numbering[i] = id
            nxyz[id, :] = fens.xyz[i, :]
            id = id + 1
        end
    end
    fens = FENodeSet(nxyz[1:(id-1), :])
    return fens, vec(new_numbering)
end

"""
    mergemeshes(
        fens1::FENodeSet{T},
        fes1::T1,
        fens2::FENodeSet{T},
        fes2::T2,
        tolerance::T,
    ) where {T, T1<:AbstractFESet, T2<:AbstractFESet}

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
function mergemeshes(
    fens1::FENodeSet{T},
    fes1::T1,
    fens2::FENodeSet{T},
    fes2::T2,
    tolerance::T,
) where {T,T1<:AbstractFESet,T2<:AbstractFESet}
    # Fuse the nodes
    # @code_warntype fusenodes(fens1, fens2, tolerance);
    fens, new_indexes_of_fens1_nodes = fusenodes(fens1, fens2, tolerance)
    # Renumber the finite elements
    newfes1 = deepcopy(fes1)
    updateconn!(newfes1, new_indexes_of_fens1_nodes)
    # Note that now the connectivity of both fes1 and fes2 point into
    # fens.
    return fens, newfes1, fes2
end

"""
    mergenmeshes(meshes::Array{Tuple{FENodeSet{T},AbstractFESet}}, tolerance::T) where {T<:Number}

Merge several meshes together.

The meshes are glued together by merging the nodes that fall within
a box of size `tolerance`. If `tolerance` is set to zero, no merging of
nodes is performed; the nodes from the meshes are simply concatenated together.

# Output
The merged node set, `fens`, and an array of finite element sets with
renumbered  connectivities are returned.
"""
function mergenmeshes(
    meshes::Array{Tuple{FENodeSet,AbstractFESet}},
    tolerance::T,
) where {T<:Number}
    outputfes = Array{AbstractFESet,1}()
    if (length(meshes)) == 1 # A single mesh, package output and return
        fens, fes = meshes[1]
        push!(outputfes, fes)
        return fens, outputfes
    end
    # Multiple meshes: process
    fens, fes = meshes[1]
    push!(outputfes, fes)
    for j = 2:length(meshes)
        fens1, fes1 = meshes[j]
        fens, new_indexes_of_fens1_nodes = fusenodes(fens1, fens, tolerance)
        updateconn!(fes1, new_indexes_of_fens1_nodes)
        push!(outputfes, fes1)
    end
    return fens, outputfes
end

"""
    mergenodes(fens::FENodeSet{T}, fes::AbstractFESet, tolerance::T) where {T<:Number}

Merge together nodes of a single node set.

Merge by gluing together nodes from a single node set located within
tolerance of each other. The nodes are glued together by merging the
nodes that fall within a box of size `tolerance`. The merged node
set, `fens`, and the finite element set, `fes`, with renumbered  connectivities
are returned.

Warning: This tends to be an expensive operation!
"""
function mergenodes(fens::FENodeSet{T}, fes::AbstractFESet, tolerance::T) where {T<:Number}
    maxnn = count(fens) + 1
    xyz1 = fens.xyz
    dim = size(xyz1, 2)
    id1 = collect(eachindex(fens))
    d = zeros(size(xyz1, 1))
    # Mark nodes from the array that are duplicated
    for i in eachindex(fens)
        if (id1[i] > 0) # This node has not yet been marked for merging
            XYZ = reshape(xyz1[i, :], 1, dim)
            copyto!(d, sum(abs.(xyz1 .- XYZ), dims = 2)) #find the distances along  coordinate directions
            minn = maxnn
            @inbounds for jx in eachindex(d)
                if d[jx] < tolerance
                    minn = min(jx, minn)
                    id1[jx] = -minn
                    id1[minn] = minn
                end
            end
        end
    end
    # Generate  merged arrays of the nodes
    xyzm = zeros(T, count(fens), dim)
    mid = 1
    for i in eachindex(fens) # and then we pick only non-duplicated fens1
        if id1[i] > 0 # this node is the master
            id1[i] = mid
            xyzm[mid, :] = xyz1[i, :]
            mid = mid + 1
        else # this node is the slave
            id1[i] = id1[-id1[i]]
        end
    end
    nnodes = mid - 1
    xyzm = xyzm[1:nnodes, :]
    # Renumber the cells
    conns = connasarray(fes)
    for i in axes(conns, 1)
        conns[i, :] = id1[conns[i, :]]
    end
    fes = fromarray!(fes, conns)

    fens = FENodeSet(xyzm[1:nnodes, :])

    return fens, fes
end

"""
    mergenodes(
        fens::FENodeSet{T},
        fes::AbstractFESet,
        tolerance::T,
        candidates::AbstractVector{IT},
    ) where {T<:Number, IT<:Integer}

Merge together  nodes of a single node set.

Similar to `mergenodes(fens, fes, tolerance)`,
but only the candidate nodes are considered for merging. This can potentially
speed up the operation by orders of magnitude.
"""
function mergenodes(
    fens::FENodeSet{T},
    fes::AbstractFESet,
    tolerance::T,
    candidates::AbstractVector{IT},
) where {T<:Number,IT<:Integer}
    maxnn = count(fens) + 1
    xyz1 = fens.xyz
    dim = size(xyz1, 2)
    id1 = collect(eachindex(fens))
    d = fill(100.0 * tolerance, size(xyz1, 1))
    # Mark nodes from the array that are duplicated
    for ic in eachindex(candidates)
        i = candidates[ic]
        if (id1[i] > 0) # This node has not yet been marked for merging
            XYZ = xyz1[i, :]
            minn = maxnn
            for kx in candidates
                d[kx] = sum(abs.(xyz1[kx, :] .- XYZ))
            end
            @inbounds for jx in candidates
                if d[jx] < tolerance
                    minn = min(jx, minn)
                    id1[jx] = -minn
                    id1[minn] = minn
                end
            end
        end
    end
    # Generate  merged arrays of the nodes
    xyzm = zeros(T, count(fens), dim)
    mid = 1
    for i in eachindex(fens) # and then we pick only non-duplicated fens1
        if id1[i] > 0 # this node is the master
            id1[i] = mid
            xyzm[mid, :] = xyz1[i, :]
            mid = mid + 1
        else # this node is the slave
            id1[i] = id1[-id1[i]]
        end
    end
    nnodes = mid - 1
    xyzm = xyzm[1:nnodes, :]
    # Renumber the cells
    conns = connasarray(fes)
    for i in eachindex(fes)
        conn = conns[i, :]
        conns[i, :] = id1[conn]
    end
    fes = fromarray!(fes, conns)

    fens = FENodeSet(xyzm[1:nnodes, :])

    return fens, fes
end

"""
    renumberconn!(fes::AbstractFESet, new_numbering::AbstractVector{IT}) where {IT<:Integer}

Renumber the nodes in the connectivity of the finite elements based on a new
numbering for the nodes.

`fes` =finite element set
`new_numbering` = new serial numbers for the nodes.  The connectivity
          should be changed as `conn[j]` --> `new_numbering[conn[j]]`

Let us say there are nodes not connected to any finite element that we would
like to remove from the mesh: here is how that would be accomplished.
```
connected = findunconnnodes(fens, fes);
fens, new_numbering = compactnodes(fens, connected);
fes = renumberconn!(fes, new_numbering);
```
Finally, check that the mesh is valid:
```julia
validate_mesh(fens, fes);
```
"""
function renumberconn!(
    fes::AbstractFESet,
    new_numbering::AbstractVector{IT},
) where {IT<:Integer}
    conn = connasarray(fes)
    for i in axes(conn, 1)
        c = conn[i, :]
        conn[i, :] = new_numbering[c]
    end
    return fromarray!(fes, conn)
end

"""
    vsmoothing(v::Matrix{T}, t::Matrix{IT}; kwargs...) where {T<:Number, IT<:Integer}

Internal routine for mesh smoothing.

Keyword options:
`method` = :taubin (default) or :laplace
`fixedv` = Boolean array, one entry per vertex: is the vertex immovable (true)
    or movable  (false)
`npass` = number of passes (default 2)
"""
function vsmoothing(v::Matrix{T}, t::Matrix{IT}; kwargs...) where {T<:Number,IT<:Integer}
    fixedv = falses(size(v, 1))
    npass = 2
    method = :taubin
    for apair in pairs(kwargs)
        sy, val = apair
        if sy == :method
            method = val
        elseif sy == :fixedv
            fixedv .= val
        elseif sy == :npass
            npass = val
        end
    end

    nv = deepcopy(v)
    # find neighbors for the given connections
    vneigh = vertexneighbors(t, size(v, 1))
    # Smoothing considering all connections through the volume
    if (method == :taubin)
        nv = smoothertaubin(v, vneigh, fixedv, npass, 0.5, -0.5)
    elseif (method == :laplace)
        nv = smootherlaplace(v, vneigh, fixedv, npass, 0.5, -0.5)
    end
    # return new vertex locations
    return nv
end

"""
    meshsmoothing(fens::FENodeSet, fes::T; options...) where {T<:AbstractFESet}

General smoothing of meshes.

# Keyword options
`method` = :taubin (default) or :laplace
`fixedv` = Boolean array, one entry per vertex: is the vertex immovable (true)
    or movable  (false)
`npass` = number of passes (default 2)

# Return
The modified  node set.
"""
function meshsmoothing(fens::FENodeSet, fes::T; options...) where {T<:AbstractFESet}
    v = deepcopy(fens.xyz)
    v = vsmoothing(v, connasarray(fes); options...)
    copyto!(fens.xyz, v)
    return fens
end

function smoothertaubin(
    vinp::Matrix{T},
    vneigh::Array{IT,1},
    fixedv::VF,
    npass,
    lambda::T,
    mu::T,
) where {T,IT,VF}
    v = deepcopy(vinp)
    nv = deepcopy(vinp)
    for I = 1:npass
        o = randperm(length(vneigh))
        damping_factor = lambda
        for k in eachindex(vneigh)
            r = o[k]
            n = vneigh[r]
            if (length(n) > 1) && (!fixedv[r])
                ln1 = (length(n) - 1)
                nv[r, :] .=
                    (1 - damping_factor) * vec(v[r, :]) +
                    damping_factor * (vec(sum(v[n, :], dims = 1)) - vec(v[r, :])) / ln1
            end
        end
        v = deepcopy(nv)
        damping_factor = mu
        for k in eachindex(vneigh)
            r = o[k]
            n = vneigh[r]
            if (length(n) > 1) && (!fixedv[r])
                ln1 = (length(n) - 1)
                nv[r, :] .=
                    (1 - damping_factor) * vec(v[r, :]) +
                    damping_factor * (vec(sum(v[n, :], dims = 1)) - vec(v[r, :])) / ln1
            end
        end
        v = deepcopy(nv)
    end
    return nv
end

function smootherlaplace(
    vinp::Matrix{T},
    vneigh::Array{IT,1},
    fixedv::VF,
    npass,
    lambda::T,
    mu::T,
) where {T,IT,VF}
    v = deepcopy(vinp)
    nv = deepcopy(vinp)
    damping_factor = lambda
    for I = 1:npass
        o = randperm(length(vneigh))
        for k in eachindex(vneigh)
            r = o[k]
            n = vneigh[r]
            if (length(n) > 1) && (!fixedv[r])
                ln1 = (length(n) - 1)
                nv[r, :] =
                    (1 - damping_factor) * vec(v[r, :]) +
                    damping_factor * (vec(sum(v[n, :], dims = 1)) - vec(v[r, :])) / ln1
            end
        end
        v = deepcopy(nv)
    end
    return nv
end

"""
    vertexneighbors(conn::Matrix{IT}, nvertices::IT) where {IT<:Integer}

Find the node neighbors in the mesh.

Return an array of integer vectors, element I holds an array of numbers of nodes
which are connected to node I (including node I).
"""
function vertexneighbors(conn::Matrix{IT}, nvertices::IT) where {IT<:Integer}
    vn = Vector{IT}[]
    sizehint!(vn, nvertices)
    for I = 1:nvertices
        push!(vn, IT[])          # preallocate
    end
    for I in axes(conn, 1)
        for r in axes(conn, 2)
            append!(vn[conn[I, r]], vec(conn[I, :]))
        end
    end
    for I in eachindex(vn)
        vn[I] = unique(vn[I])
    end
    return vn
end

"""
    mirrormesh(
        fens::FENodeSet,
        fes::ET,
        Normal::Vector{T},
        Point::Vector{T};
        kwargs...,
    ) where {ET<:AbstractFESet, T<:Number}

Mirror a mesh in a plane given by its normal and one point.

# Keyword arguments
- `renumb` = node renumbering function for the mirrored element

Warning: The code to relies on the numbering of the finite elements: to reverse
the orientation of the mirrored finite elements, the connectivity is listed in
reverse order.   If the mirrored finite elements do not follow this rule (for
instance hexahedra or quadrilaterals), their areas/volumes will come out
negative. In such a case the renumbering function of the connectivity needs to
be supplied.

For instance: H8 elements require the renumbering function to be supplied as
```
renumb = (c) -> c[[1, 4, 3, 2, 5, 8, 7, 6]]
```
"""
function mirrormesh(
    fens::FENodeSet,
    fes::ET,
    Normal::Vector{T},
    Point::Vector{T};
    kwargs...,
) where {ET<:AbstractFESet,T<:Number}
    # Default renumbering function.
    # Simply switch the order of nodes.  Works for simplexes...
    renumb(conn) = conn[end:-1:1]
    for apair in pairs(kwargs)
        sy, val = apair
        if sy == :renumb
            renumb = val
        end
    end
    # Make sure we're using a unit normal
    Normal = Normal / norm(Normal)
    Normal = vec(Normal)
    # The point needs to be a row  matrix
    Point = vec(Point)

    fens1 = deepcopy(fens) # the mirrored mesh nodes
    for i in eachindex(fens1)
        a = fens1.xyz[i, :]
        d = dot(vec(a - Point), Normal)
        fens1.xyz[i, :] = a - 2 * d * Normal
    end
    # Reconnect the cells
    fes1 = deepcopy(fes)
    conn = connasarray(fes1)
    for i in axes(conn, 1)
        conn[i, :] = renumb(conn[i, :])
    end
    return fens1, fromarray!(fes1, conn)
end

function _nodepartitioning3(xyz, nincluded::Vector{Bool}, npartitions::Int = 2)
    function inertialcutpartitioning!(partitioning, parts, X)
        nspdim = 3
        StaticMoments = fill(zero(eltype(xyz)), nspdim, length(parts))
        npart = fill(0, length(parts))
        for spdim = 1:nspdim
            @inbounds for j in axes(X, 1)
                if nincluded[j] # Is the node to be included in the partitioning?
                    jp = partitioning[j]
                    StaticMoments[spdim, jp] += X[j, spdim]
                    npart[jp] += 1 # count the nodes in the current partition
                end
            end
        end
        CG = fill(zero(eltype(xyz)), nspdim, length(parts))
        for p in parts
            npart[p] = Int(npart[p] / nspdim)
            CG[:, p] = StaticMoments[:, p] / npart[p] # center of gravity of each partition
        end
        MatrixMomentOfInertia = fill(zero(eltype(xyz)), nspdim, nspdim, length(parts))
        @inbounds for j in axes(X, 1)
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
        for p in parts
            MatrixMomentOfInertia[2, 1, p] = MatrixMomentOfInertia[1, 2, p]
            MatrixMomentOfInertia[3, 1, p] = MatrixMomentOfInertia[3, 1, p]
            MatrixMomentOfInertia[3, 2, p] = MatrixMomentOfInertia[3, 2, p]
        end
        longdir = fill(zero(eltype(xyz)), nspdim, length(parts))
        for p in parts
            F = eigen(MatrixMomentOfInertia[:, :, p])
            six = sortperm(F.values)
            longdir[:, p] = F.vectors[:, six[1]]
        end
        toggle = fill(one(eltype(xyz)), length(parts))
        @inbounds for j in axes(X, 1)
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

    nlevels = Int(round(ceil(log(npartitions) / log(2))))
    partitioning = fill(1, size(xyz, 1))  # start with nodes assigned to partition 1
    for level = 0:1:(nlevels-1)
        inertialcutpartitioning!(partitioning, collect(1:(2^level)), xyz)
    end
    return partitioning
end

function _nodepartitioning2(xy, nincluded::Vector{Bool}, npartitions::Int = 2)
    function inertialcutpartitioning!(partitions, parts, X)
        nspdim = 2
        StaticMoments = fill(zero(eltype(xy)), nspdim, length(parts))
        npart = fill(0, length(parts))
        for spdim = 1:nspdim
            @inbounds for j in axes(X, 1)
                if nincluded[j] # Is the node to be included in the partitioning?
                    jp = partitions[j]
                    StaticMoments[spdim, jp] += X[j, spdim]
                    npart[jp] += 1 # count the nodes in the current partition
                end
            end
        end
        CG = fill(zero(eltype(xy)), nspdim, length(parts))
        for p in parts
            npart[p] = Int(npart[p] / nspdim)
            CG[:, p] = StaticMoments[:, p] / npart[p] # center of gravity of each partition
        end
        MatrixMomentOfInertia = fill(zero(eltype(xy)), nspdim, nspdim, length(parts))
        @inbounds for j in axes(X, 1)
            if nincluded[j] # Is the node to be included in the partitioning?
                jp = partitions[j]
                xj, yj = X[j, 1] - CG[1, jp], X[j, 2] - CG[2, jp]
                MatrixMomentOfInertia[1, 1, jp] += yj^2
                MatrixMomentOfInertia[2, 2, jp] += xj^2
                MatrixMomentOfInertia[1, 2, jp] -= xj * yj
            end
        end
        for p in parts
            MatrixMomentOfInertia[2, 1, p] = MatrixMomentOfInertia[1, 2, p]
        end
        longdir = fill(zero(eltype(xy)), nspdim, length(parts))
        for p in parts
            F = eigen(MatrixMomentOfInertia[:, :, p])
            six = sortperm(F.values)
            longdir[:, p] = F.vectors[:, six[1]]
        end
        toggle = fill(one(eltype(xy)), length(parts))
        @inbounds for j in axes(X, 1)
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

    nlevels = Int(round(ceil(log(npartitions) / log(2))))
    partitions = fill(1, size(xy, 1))  # start with nodes assigned to partition 1
    for level = 0:1:(nlevels-1)
        inertialcutpartitioning!(partitions, collect(1:(2^level)), xy)
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
  File =  "partition-nodes-\$(gp).vtk"
  vtkexportmesh(File, fens, FESetP1(reshape(groupnodes, length(groupnodes), 1)))
end
File =  "partition-mesh.vtk"
vtkexportmesh(File, fens, fes)
@async run(`"paraview.exe" \$File`)
```
"""
function nodepartitioning(fens::FENodeSet, nincluded::Vector{Bool}, npartitions::Int = 2)
    (npartitions >= 2) || error("Number of partitions must be >= 2")
    if size(fens.xyz, 2) == 3
        return _nodepartitioning3(fens.xyz, nincluded, npartitions)
    elseif size(fens.xyz, 2) == 2
        return _nodepartitioning2(fens.xyz, nincluded, npartitions)
    else
        @warn "Not implemented for 1D"
    end
end

"""
    nodepartitioning(fens::FENodeSet, npartitions)

Compute the inertial-cut partitioning of the nodes.

`npartitions` = number of partitions, but note that the actual number of
partitions will be a power of two.

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
    (npartitions >= 2) || error("Number of partitions must be >= 2")
    nincluded = fill(true, count(fens)) # The default is all nodes are included in the partitioning.
    if size(fens.xyz, 2) == 3
        return _nodepartitioning3(fens.xyz, nincluded, npartitions)
    elseif size(fens.xyz, 2) == 2
        return _nodepartitioning2(fens.xyz, nincluded, npartitions)
    else
        @warn "Not implemented for 1D"
    end
end

"""
    pointpartitioning(xyz, npartitions::Int = 2)

Compute the inertial-cut partitioning of a set of points.
"""
function pointpartitioning(xyz, npartitions::Int = 2)
    (npartitions >= 2) || error("Number of partitions must be >= 2")
    nincluded = fill(true, size(xyz, 1)) # The default is all nodes are included in the partitioning.
    if size(xyz, 2) == 3
        return _nodepartitioning3(xyz, nincluded, npartitions)
    elseif size(xyz, 2) == 2
        return _nodepartitioning2(xyz, nincluded, npartitions)
    else
        @warn "Not implemented for 1D"
    end
end

"""
    nodepartitioning(fens::FENodeSet, fesarr, npartitions::Vector{Int})

Compute the inertial-cut partitioning of the nodes.

`fesarr` = array of finite element sets that represent regions
`npartitions` = array of the number of partitions in each region. However
note that the actual number of partitions will be a power of two.

The partitioning itself is carried out by `nodepartitioning()` with
a list of nodes to be included in the partitioning. For each region I
the nodes included in the partitioning are those connected to
the elements of that region, but not to elements that belong to
any of the previous regions, 1…I-1.
"""
function nodepartitioning(fens::FENodeSet, fesarr, npartitions::Vector{Int})
    (length(fesarr) == length(npartitions)) || error("Number of regions and number of partitions must match")
    # Partitioning of all the nodes
    partitioning = fill(0, count(fens))
    # Find the partitioning of the nodes in FESet 1
    nincludedp = fill(false, count(fens))
    for i in connectednodes(fesarr[1]) # For nodes connected by region 1
        nincludedp[i] = true
    end
    partitioning1 = nodepartitioning(fens, nincludedp, npartitions[1])
    totnpartitions = maximum(partitioning1)
    # Transfer the partitioning of region 1 into the overall partitioning
    partitioning[nincludedp] = partitioning1[nincludedp]
    for i = 2:length(npartitions)
        # Find the partitioning of the nodes in FESet i, but not in the preceding sets
        nincluded = fill(false, count(fens))
        for j in connectednodes(fesarr[i]) # For nodes connected by region i
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

"""
    distortblock(ofens::FENodeSet{T}, xdispmul::T, ydispmul::T) where {T<:Number}

Distort a block mesh by shifting around the nodes. The goal is to
distort the horizontal and vertical mesh lines into slanted lines.
"""
function distortblock(ofens::FENodeSet{T}, xdispmul::T, ydispmul::T) where {T<:Number}
    Lx = maximum(ofens.xyz[:, 1])
    Ly = maximum(ofens.xyz[:, 2])
    xic(x) = (2 * x - Lx) / Lx
    etac(y) = (2 * y - Ly) / Ly

    fens = deepcopy(ofens)

    for k in eachindex(fens)
        x, y = fens.xyz[k, :]
        xi, eta = xic(x), etac(y)
        u = xdispmul * (1 + xi) * (1 - xi) * eta
        v = ydispmul * xi * (1 + eta) * (1 - eta)

        fens.xyz[k, 1] = x + u
        fens.xyz[k, 2] = y + v
    end

    return fens
end

"""
    distortblock(
        B::F,
        Length::T,
        Width::T,
        nL::IT,
        nW::IT,
        xdispmul::T,
        ydispmul::T,
    ) where {F<:Function, T<:Number, IT<:Integer}

Distort a block mesh by shifting around the nodes. The goal is to distort the
horizontal and vertical mesh lines into slanted lines. This is useful when
testing finite elements where special directions must be avoided.
"""
function distortblock(
    B::F,
    Length::T,
    Width::T,
    nL::IT,
    nW::IT,
    xdispmul::T,
    ydispmul::T,
) where {F<:Function,T<:Number,IT<:Integer}
    (1.0 >= abs(xdispmul) >= 0.0) || error("xdispmul must be between 0 and 1")
    (1.0 >= abs(ydispmul) >= 0.0) || error("ydispmul must be between 0 and 1")
    fens, fes = B(1.0, 1.0, nL, nW)
    nfens = deepcopy(fens)
    if xdispmul != 0.0 || ydispmul != 0.0
        nfens = distortblock(fens, xdispmul, ydispmul)
    end
    nfens.xyz[:, 1] .*= Length
    nfens.xyz[:, 2] .*= Width
    return nfens, fes
end

function __all_behind(xyz, centroid, c, n, conn, tol, mask)
    v = deepcopy(c)
    #     Check the centroid first
    @. v = centroid - c
    vn = norm(v)
    if (vn > 0)
        d = dot(n, v) / vn
        if (d > tol)
            return false
        end
    end
    # Now check all the vertices
    # except for the nodes which are  connected by the tested cell
    mask .= true
    mask[conn] .= false
    for k in eachindex(mask)
        if mask[k]
            v .= vec(view(xyz, k, :)) - c
            vn = norm(v)
            if (vn > 0)
                d = dot(n, v) / vn
                if (d > tol)
                    return false
                end
            end
        end
    end
    return true
end

"""
    outer_surface_of_solid(fens, bdry_fes)

Find the finite elements that form the outer boundary surface.

!!! note:

The code will currently not work correctly for 2D axially symmetric geometries.

# Return

Set of boundary finite elements that enclose the solid. Now cavities 
are included.
"""
function outer_surface_of_solid(fens::FENodeSet, bdry_fes::ET) where {ET<:AbstractFESet}
    sn = SurfaceNormal(size(fens.xyz, 2))
    parametric_centroid = centroidparametric(bdry_fes)
    N = bfun(bdry_fes, parametric_centroid)
    gradNpar = bfundpar(bdry_fes, parametric_centroid)
    # This is the centroid of the cloud of nodes
    centroid = vec(mean(fens.xyz, dims = 1))
    angtol = 0.001

    # # for axially symmetric geometry place the centroid on the axis
    # if (get(bdry_fes,'axisymm'))
    #     centroid(1)=0;
    # end

    # If we have 2-D axial symmetry: eliminate  boundary cells that are on the
    # axis of symmetry
    # onaxisl  = [];
    # if (get(bdry_fes,'axisymm'))
    #     onaxisl = fe_select(fens,bdry_fes,struct ('box',[0,0,-inf,inf],'inflate',norm(max(xyz)-min(xyz))/10000));
    #     # Prune the  list of cells to include by flooding below
    #     bdry_fes=subset(bdry_fes,setdiff(eachindex(bdry_fes),onaxisl));;
    # end

    # Now we will try to find  one  boundary cell  so that  all the nodes lie in
    # the half space  behind it (that is in the other half space from the one
    # into which  its  normal is pointing), and the half space also includes
    # the centroid.
    conns = connasarray(bdry_fes)
    mask = fill(true, count(fens))
    start = 0
    for i in eachindex(bdry_fes)
        c = N' * view(fens.xyz, conns[i, :], :)
        J = view(fens.xyz, conns[i, :], :)' * gradNpar
        n = updatenormal!(sn, c, J, i, 0)
        n ./= norm(n)
        if __all_behind(fens.xyz, centroid, vec(c), n, conns[i, :], angtol, mask)
            start = i
            break
        end
    end

    # Could not find a single  surface cell  so that all nodes are behind it:
    # failure.
    if start == 0
        return nothing
    end

    startfen = conns[start, 1]

    # Select all cells connected together
    osfesl = selectelem(fens, bdry_fes, flood = true, startnode = startfen)

    return subset(bdry_fes, osfesl)
end

"""
    reordermesh(fens, fes, ordering)    

Reorder mesh (reshuffle nodes, renumber connectivities correspondingly).

# Arguments
- `fens`: The set of mesh nodes.
- `fes`: The set of elements.
- `ordering`: The desired ordering of the nodes and elements.

# Returns
The reordered mesh nodes and elements.

The ordering may come from Reverse Cuthill-McKey (package SymRCM).
"""
function reordermesh(fens, fes, ordering)
    iordering = collect(1:length(ordering))
    iordering[ordering] = iordering
    return FENodeSet(fens.xyz[ordering, :]), renumberconn!(fes, iordering)
end

"""
    element_coloring(fes, n2e, ellist::Vector{IT} = IT[]) where {IT<:Integer}

Find coloring of the elements such that no two elements of the same color share
a node.

Returns the colors of the elements (color here means an integer), and a list
of the unique colors (numbers).

`ellist` = list of elements to be assigned colors; other element colors may be looked at
"""
function element_coloring(fes, n2e, ellist::Vector{IT} = Int[]) where {IT<:Integer}
    element_colors = fill(zero(eltype(fes.conn[1])), count(fes))
    unique_colors = eltype(element_colors)[1]
    color_counts = eltype(element_colors)[0]
    if isempty(ellist) 
        ellistit = 1:count(fes)
    else
        ellistit = ellist
    end
    return element_coloring!(element_colors, unique_colors, color_counts, fes, n2e, ellistit)
end

function __find_minimal_count(color_used, color_counts)
    c = 0
    mincount = typemax(eltype(color_counts))
    @inbounds for k in eachindex(color_used)
        if color_used[k] == 0 && mincount > color_counts[k]
            mincount = color_counts[k]
            c = k
        end
    end
    return c
end

function __find_first_zero(color_used)
    @inbounds for k in eachindex(color_used)
        if color_used[k] == 0 
            return k
        end
    end
    return 0
end

"""
    element_coloring!(element_colors, unique_colors, color_counts, fes, n2e, ellist)

Find coloring of the elements such that no two elements of the same color share
a node.

Compute element coloring, supplying the current element colors and the list of
elements to be assigned colors.
"""
function element_coloring!(element_colors, unique_colors, color_counts, fes, n2e, ellist_iterator)
    __color_used = fill(zero(eltype(element_colors)), length(unique_colors))
    for e in ellist_iterator
        if element_colors[e] == 0
            __color_used .= 0
            @inbounds for n in fes.conn[e]
                m = n2e.map[n]
                for j in eachindex(m)
                    oe = m[j]
                    c = element_colors[oe]
                    if c > 0
                        __color_used[c] += 1
                    end
                end
            end
            if sum(__color_used) == 0
                c = argmin(color_counts)
                element_colors[e] = unique_colors[c]
                color_counts[c] += 1
            else
                if __find_first_zero(__color_used) == 0
                    added = maximum(unique_colors) + 1
                    push!(unique_colors, added)
                    push!(color_counts, 0)
                    push!(__color_used, 0)
                    c = length(color_counts)
                    element_colors[e] = unique_colors[c]
                    color_counts[c] += 1
                else
                    c = __find_minimal_count(__color_used, color_counts)
                    element_colors[e] = unique_colors[c]
                    color_counts[c] += 1
                end
            end
        end
    end
    return element_colors, unique_colors, color_counts
end

"""
    validate_mesh(fens, fes)

Validate the given mesh by checking if it satisfies certain sanity criteria.

# Arguments
- `fens`: The finite element nodes of the mesh.
- `fes`: The finite elements of the mesh.

Validate finite element mesh.
    
A finite element mesh given by the node set and the finite element set is
validated by checking the sanity of the numbering:
- the node numbers need to be positive and in serial order
- the fe connectivity needs to refer to valid nodes
- the finite element nodes need to be connected to at least one finite element

An error is reported as soon as it is detected.

# Returns
A boolean indicating whether the mesh is valid or not.
"""
function validate_mesh(fens, fes)
    totnfens = count(fens)
    for i in eachindex(fes)
        if (max(fes.conn[i]...) > totnfens)
            error("Wrong connectivity (refers to nonexistent node): $(fes.conn[i])")
        end
        if (min(fes.conn[i]...) < 1)
            error("Wrong connectivity (refers to nonexistent node): $(fes.conn[i])")
        end
        if (length(unique(fes.conn[i])) != length(fes.conn[i]))
            error("Wrong connectivity (multiply referenced node): $(fes.conn[i])")
        end
    end
    connected = findunconnnodes(fens, fes)
    if (any(connected == 0))
        error("Unconnected nodes present: $(sum(connected==0)) total")
    end
    return true
end


end # module
