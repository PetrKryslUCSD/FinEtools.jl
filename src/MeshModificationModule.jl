module MeshModificationModule

using FinEtools.FTypesModule
using FinEtools.FENodeSetModule
using FinEtools.FESetModule
using Base.Sort
using Base.Order

"""
    meshboundary{T<:FESet}(fes::T)

Extract the boundary finite elements from a mesh.

Extract the finite elements of manifold dimension (n-1) from the
supplied list of finite elements of manifold dimension (n).
"""
function meshboundary{T<:FESet}(fes::T)
    # Form all hyperfaces, non-duplicates are boundary cells
    hypf= FESetModule.boundaryconn(fes);    # get the connectivity of the boundary elements
    bdryconn = myunique2(hypf);
    make = FESetModule.boundaryfe(fes);     # get the function that can make a boundary element

    return make(bdryconn);
end
export meshboundary

function mysortrows(A::FIntMat)
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

function mysortdim2!(A::FIntMat)
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

function  myunique2(A::FIntMat) # speeded up; now the bottleneck is mysortrows
    #println("size(A)=$(size(A))")
    maxA=maximum(A[:])::FInt
    sA=deepcopy(A)
    #@time
    sA=mysortdim2!(sA)::FIntMat;#this is fast
    #@time sA=sort(A,2,alg=QuickSort)::FIntMat;#this is slow
    sA= [sA (1:size(A,1))+maxA]::FIntMat
    #@time
    sA =mysortrows(sA); # this now takes the majority of time, but much less than the function below
    #@time sA  = sortrows(sA,alg=QuickSort);;#this is slow
    rix=sA[:,end]; rix=rix[:]-maxA;
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

### This code is correct, but very slow.
function  myunique1(A::FIntMat)
    maxA=maximum(A[:])
    sA=sort(A,2);# most time spent here
    sA= [sA (1:size(A,1))+maxA]
    sA  = sortrows(sA);;#and here
    rix=sA[:,end]; rix=rix[:]-maxA;
    sA=sA[:,1:end-1];
    d=(sA[1:end-1,:].!=sA[2:end,:]); # element-wise comparison!
    ad=map((x) -> (x?1:0),[true; any(d,2)]);
    iu =map((x) -> (x>1? true: false),(ad + [ad[2:end];1]));
    Out =A[rix[iu[:]],:];
    return Out
end

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
k=new_indexes_of_fens1_nodes[j] is the node in the node set fens which
used to be node `j` in node set `fens1`.
The finite element set connectivity that used to refer to `fens1`
needs to be updated to refer to the same nodes in  the set `fens` as
     `updateconn!(fes, new_indexes_of_fens1_nodes);`
"""
function fusenodes(fens1::FENodeSet, fens2::FENodeSet, tolerance:: FFlt)
  xyz1::FFltMat = deepcopy(fens1.xyz);
  id1::FIntVec = collect(1:size(xyz1,1));
  dim =size(xyz1,2);
  xyz2::FFltMat = deepcopy(fens2.xyz);
  id2::FIntVec = collect(1:size(xyz2,1));
  n1::FFlt= 0.0
  # Mark nodes from the first array that are duplicated in the second
  if (tolerance>0.0) # should we attempt to merge nodes?
    for i=1:size(xyz1,1)
      for rx=1:size(xyz2,1)
        n1= 0.0
        for cx=1:size(xyz2,2)
          n1=n1+abs(xyz2[rx,cx]-xyz1[i,cx]);
        end
        if (n1<tolerance)
          id1[i] =-rx; break;
        end
      end
    end
  end
  # Generate  fused arrays of the nodes
  xyzm = zeros(FFlt,size(xyz1,1)+size(xyz2,1),dim);
  for rx=1:size(xyz2,1)
    xyzm[rx,:]=xyz2[rx,:];
  end
  idm = zeros(FInt,size(xyz1,1)+size(xyz2,1));
  for rx=1:size(xyz2,1)
    idm[rx]=rx;
  end
  mid=size(xyz2,1)+1;
  for i=1:size(xyz1,1) # and then we pick only non-duplicated fens1
    if id1[i]>0
      id1[i]=mid;
      idm[mid]=mid;
      for cx=1:size(xyz1,2)
        xyzm[mid,cx]=xyz1[i,cx];
      end
      mid=mid+1;
    else
      id1[i]=id2[-id1[i]];
    end
  end
  nnodes =mid-1;
  xyzm =xyzm[1:nnodes,:];

  # Create the fused Node set
  fens = FENodeSet(xyzm);
  # The Node set 1 numbering will change
  new_indexes_of_fens1_nodes = id1[:];
  # The node set 2 numbering stays the same
  return fens, new_indexes_of_fens1_nodes
end
export fusenodes

"""
    compactnodes(fens::FENodeSetModule.FENodeSet, connected::Vector{Bool})

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

connected = findunconnnodes(fens, fes);
fens, new_numbering =compactnodes(fens, connected);
fes = renumberconn!(fes, new_numbering);

Finally, check that the mesh is valid:
validate_mesh(fens, fes);
"""
function compactnodes(fens::FENodeSetModule.FENodeSet, connected::BitArray{1})
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
  fens = FENodeSetModule.FENodeSet(nxyz[1:id-1,:]);
  return fens, vec(new_numbering)
end
export compactnodes

"""
    mergemeshes(fens1::FENodeSet, fes1::T1,
      fens2::FENodeSet, fes2::T2, tolerance::FFlt) where {T1<:FESet,T2<:FESet}

Merge together two meshes.

Merge two meshes together by gluing together nodes within tolerance. The
two meshes, fens1, fes1, and fens2, fes2, are glued together by merging
the nodes that fall within a box of size "tolerance". If tolerance is set
to zero, no merging of nodes is performed; the two meshes are simply
concatenated together.

The merged node set, fens, and the two arrays of finite elements with
renumbered  connectivities are returned.

Important notes: On entry into this function the connectivity of fes1
point into fens1 and the connectivity of fes2 point into fens2. After
this function returns the connectivity of both fes1 and fes2 point into
fens. The order of the nodes of the node set fens1 in the resulting set
fens will have changed, whereas the order of the nodes of the node set
fens2 is are guaranteed to be the same. Therefore, the connectivity of
fes2 will in fact remain the same.
"""
function mergemeshes(fens1::FENodeSet, fes1::T1,
  fens2::FENodeSet, fes2::T2, tolerance::FFlt) where {T1<:FESet,T2<:FESet}
  # Fuse the nodes
  fens, new_indexes_of_fens1_nodes = fusenodes(fens1, fens2, tolerance);
  # Renumber the finite elements
  newfes1 = deepcopy(fes1)
  updateconn!(newfes1, new_indexes_of_fens1_nodes);
  # Note that now the connectivity of both fes1 and fes2 point into
  # fens.
  return fens, newfes1, fes2
end
export mergemeshes

function mergenmeshes!(fensa, fesa, tolerance::FFlt)
#     Merge several meshes together.
# %
# function [fens,fesa] = merge_n_meshes(fensa, fesa, tolerance)
# %
# Merge several meshes together either by simple concatenation of nodes or by
# gluing together nodes within tolerance.
# %
# Inputs:
# fensa= cell array of node sets, one for each mesh;
# fesa= cell array of finite element sets, one for each mesh;
# tolerance= Geometric tolerance, maybe supplied as zero (>=0).
# %
# The meshes are glued together by
# merging the nodes that fall within a box of size "tolerance". If tolerance
# is set to zero, no merging of nodes is performed; the nodes from the meshes are
# simply concatenated together.
# %
# The merged node set, fens, and the cell array of finite element sets with
# renumbered  connectivities are returned.
# %
# Outputs:
# fens= merged node set,
# fesa= cell array of finite element sets updated to use the merged node set.
# %
# %
# See also: merge_meshes


    if (length(fensa))!=(length(fesa))
        error("(length(fensa))!=(length(fesa))");
    end
    if (length(fensa))==1
        fens=fensa[1];
        return fens,fesa        # There is nothing to be done: this is a single mesh
    end
    fens=fensa[1];
    for j=2:length(fesa)
        fens,new_indexes_of_fens1_nodes = fusenodes(fensa[j], fens, tolerance);
        FESetModule.updateconn!(fesa[j],new_indexes_of_fens1_nodes);
    end
    return fens,fesa
end
export mergenmeshes!

"""
    mergenodes(fens::FENodeSet, fes::FESet, tolerance::FFlt)

Merge together  nodes of a single node set.

Merge by gluing together nodes from a single node set located within
tolerance of each other. The nodes are glued together by merging the
nodes that fall within a box of size `tolerance`. The merged node
set, fens, and the finite element set with renumbered  connectivities
are returned.
"""
function mergenodes(fens::FENodeSet, fes::FESet, tolerance::FFlt)
  xyz1 = fens.xyz;
  dim  = size(xyz1,2);
  id1 = collect(1:count(fens));
  c1 = ones(size(xyz1,1),1);
  xyzd = zeros(size(xyz1));
  d = zeros(size(xyz1,1));
  m = trues(size(xyz1,1));
  # Mark nodes from the array that are duplicated
  for i = 1:count(fens)
    if (id1[i]>0) # This node has not yet been marked for merging
      XYZ = reshape(xyz1[i,:], 1, dim);
      xyzd[:,:] = abs.(xyz1-c1*XYZ); #find the distances along  coordinate directions
      d = sum(xyzd,2);
      map!((x)->x<tolerance, m, d);
      jx = find(m);
      if (!isempty(jx))
        minn = minimum(jx);
        id1[jx] = -minn;
        id1[minn] = minn;
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
  conns = fes.conn;
  for i = 1:FESetModule.count(fes)
    conn = conns[i,:];
    conns[i,:] = id1[conn];
  end
  fes.conn = deepcopy(conns);

  fens = FENodeSet(xyzm[1:nnodes,:]);

  return fens,fes
end
export mergenodes

"""
    renumberconn!(fes::FESetModule.FESet, new_numbering::FIntVec)

Renumber the nodes in the connectivity of the finite elements based on a new
numbering for the nodes.

fes =finite element set
new_numbering = new serial numbers for the nodes.  The connectivity
          should be changed as conn[j] --> new_numbering(conn[j])

Let us say there are nodes not connected to any finite element that you would
like to remove from the mesh: here is how that would be accomplished.
%
connected = findunconnnodes(fens, fes);
fens, new_numbering =compactfens(fens, connected);
fes = renumberconn!(fes, new_numbering);
%
Finally, check that the mesh is valid:
validate_mesh(fens, fes);
"""
function renumberconn!(fes::FESetModule.FESet, new_numbering::FIntVec)
  for i=1:size(fes.conn,1)
    c = fes.conn[i,:];
    fes.conn[i,:] = new_numbering[c];
  end
  return fes
end
export renumberconn!

function vsmoothing(v::FFltMat,t::FIntMat;options...)
# General smoothing of meshes.
#
# function [t,v] =smoothing(t,v,options)
#
# Fields of the structure options, all are optional:
# method='laplace' or 'taubin' (Default is 'taubin'.)
# f=boundary faces (optional)
# bv=boundary vertices (optional)
# bv_from_f=compute boundary of vertices from the boundary faces, true or
# false.  Tetrahedra and hexahedra are supported.
    # npass=how many passes of smoothing? default is 2.


    iv=deepcopy(v);

    fixedv=falses(size(v,1))
    npass =2;
    method =:taubin;
    for arg in options
        sy, val = arg
        if sy==:method
            method=val
        elseif sy==:fixedv
            fixedv=val
        elseif sy==:npass
            npass=val
        end
    end

    # find neighbors for the given connections
    vneigh =  vertex_neighbors(t,size(v,1));
    # Smoothing considering all connections through the volume
    if (method==:taubin)
        v =  taubin_smoother(v,vneigh,fixedv,npass,0.5,-0.5);
    elseif (method==:laplace)
        v =  laplace_smoother(v,vneigh,fixedv,npass,0.5,-0.5);
    end
    # return new vertex locations
    return v
end
export vsmoothing


function meshsmoothing{T<:FESet}(fens::FENodeSet,fes::T;options...)
# General smoothing of meshes.
#
# function [t,v] =smoothing(t,v,options)
#
# Fields of the structure options, all are optional:
# method='laplace' or 'taubin' (Default is 'taubin'.)
# f=boundary faces (optional)
# bv=boundary vertices (optional)
# bv_from_f=compute boundary of vertices from the boundary faces, true or
# false.  Tetrahedra and hexahedra are supported.
    # npass=how many passes of smoothing? default is 2.
    nnodes=deepcopy(fens)
    v=vsmoothing(nnodes.xyz,fes.conn;options...)
    nnodes.xyz= deepcopy(v)
    return nnodes,fes
end
export meshsmoothing

function  taubin_smoother(vinp::FFltMat,vneigh::Array{Array{Int,1},1},fixedv::BitArray{1},npass::FInt,lambda::FFlt,mu::FFlt)
    v=deepcopy(vinp);
    nv=deepcopy(vinp);
    for I= 1:npass
        o=randperm(length(vneigh));
        damping_factor=lambda;
        for k= 1:length(vneigh)
            r=o[k];
            n=vneigh[r];
            if (length(n)>1) && (!fixedv[r])
                nv[r,:]=(1-damping_factor)*v[r,:]+ damping_factor*(sum(v[n,:],1)-v[r,:])/(length(n)-1);
            end
        end
        v=deepcopy(nv);
        damping_factor=mu;
        for k= 1:length(vneigh)
            r=o[k];
            n=vneigh[r];
            if (length(n)>1) && (!fixedv[r])
                nv[r,:]=(1-damping_factor)*v[r,:]+ damping_factor*(sum(v[n,:],1)-v[r,:])/(length(n)-1);
            end
        end
        v=deepcopy(nv);
    end
    return nv
end

function   laplace_smoother(vinp::FFltMat,vneigh::Array{Array{Int,1},1},fixedv::BitArray{1},npass::FInt,lambda::FFlt,mu::FFlt)
    v=deepcopy(vinp);
    nv=deepcopy(vinp);
    damping_factor=lambda;
    for I= 1:npass
        o=randperm(length(vneigh));
        for k= 1:length(vneigh)
            r=o[k];
            n=vneigh[r];
            if (length(n)>1) && (!fixedv[r])
                nv[r,:]=(1-damping_factor)*v[r,:]+ damping_factor*(sum(v[n,:],1)-v[r,:])/(length(n)-1);
            end
        end
        v=deepcopy(nv);
    end
    return nv
end

function vertex_neighbors(conn::FIntMat,nvertices::FInt)
# Find the node neighbors in the mesh.
# %
# function vn =  vertex_neighbors(vn,f,v)
# %
# vn= cell array, element I holds an array of numbers of nodes
#     which are connected to node I (including node I).  When this array is
#     supplied as input the information from the current call is added to
#     the array vn; otherwise (when vn is empty on input) the array is created
#     and returned.
# f= connectivity of the mesh, one row per element
# v= locations of the nodes, three columns, one row per node
    vn=Array(FIntVec,nvertices)
    for I= 1:length(vn)
        vn[I]=FInt[];          # preallocate
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
      Point::FFltVec; args...) where {T<:FESet}

Mirror a 2-D mesh in a plane given by its normal and one point.

Warning: The code to relies on the numbering of the cells: to reverse
the orientation of the mirrored cells, the connectivity is listed in
reverse order.   If the mirrored cells do not follow this rule (for instance
hexahedra for quadrilaterals), their areas/volumes will
come out negative. In such a case the renumbering function
of the connectivity needs to be supplied.

For instance: H8 elements require  the renumbering function to be supplied as
fens1,gcells1 = mirror_mesh(fens, gcells,...
          [-1,0,0], [0,0,0], @(c)c([1, 4, 3, 2, 5, 8, 7, 6]));
"""
function mirrormesh(fens::FENodeSet, fes::T, Normal::FFltVec,
  Point::FFltVec; args...) where {T<:FESet}
  # Treat optional arguments.
  # Simply switch the order of nodes.  Works for simplexes...
  renumb(conn) = conn[end:-1:1];
  for arg in args
    sy, val = arg
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
  for i=1:size(fes1.conn,1)
    fes1.conn[i,:]=renumb(fes1.conn[i,:]);
  end
  return fens1,fes1
end
export mirrormesh


end
