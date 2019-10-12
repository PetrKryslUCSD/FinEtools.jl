"""
    MeshSelectionModule

Module for  selection of mesh entities.
"""
module MeshSelectionModule

using ..FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import ..FESetModule: AbstractFESet, bfundpar, nodesperelem, manifdim, connasarray
import ..FENodeSetModule: FENodeSet, spacedim
import ..BoxModule: inflatebox!, initbox!, updatebox!, boxesoverlap
import ..FENodeToFEMapModule: FENodeToFEMap
import LinearAlgebra: norm, dot, cross
import Statistics: mean

"""
    connectednodes(fes::AbstractFESet)

Extract the node numbers of the nodes connected by given finite elements.

Extract the list of unique node numbers for the nodes that are connected
by the finite element set `fes`. Note that it is assumed that all the
FEs are of the same type (the same number of connected nodes by each
cell).
"""
function connectednodes(fes::AbstractFESet)
    return unique(connasarray(fes)[:]);
end

"""
    connectedelems(fes::AbstractFESet, node_list::FIntVec)

Extract the list of numbers for the fes  that are connected to given
nodes.
"""
function connectedelems(fes::AbstractFESet, node_list::FIntVec, nmax::FInt)
    f2fm = FENodeToFEMap(fes.conn, nmax)
    cg = zeros(FInt, count(fes)); # No elements are part of the group to begin with
    for j = node_list # for all nodes in the list
        for i = 1:length(f2fm.map[j])
            cg[f2fm.map[j][i]] = 1   # Mark element as being part of the group
        end
    end
    return findall(v -> (v != 0), cg)
end

"""
    selectnode(fens::FENodeSet; kwargs...)

Select nodes using some criterion.

# Arguments
- `v` = array of locations, one location per row
- `kwargs` = pairs of keyword argument/value

## Selection criteria

### box
```
nLx = vselect(fens.xyz, box = [0.0 Lx  0.0 0.0 0.0 0.0], inflate = Lx/1.0e5)
```

The keyword 'inflate' may be used to increase or decrease the extent of
the box (or the distance) to make sure some nodes which would be on the
boundary are either excluded or included.

### distance
```
list = selectnode(fens.xyz, distance=1.0+0.1/2^nref, from=[0. 0.],
        inflate=tolerance);
```

### plane
```
candidates = selectnode(fens, plane = [0.0 0.0 1.0 0.0], thickness = h/1000)
```
The keyword `plane` defines the plane by its normal (the first two or
three numbers) and its distance from the origin (the last number). Nodes
are selected they lie on the plane,  or near the plane within the
distance `thickness` from the plane. The normal is assumed to be of unit
length, if it isn't apply as such, it will be normalized internally.

### nearestto
Find the node nearest to the location given.
```
nh = selectnode(fens, nearestto = [R+Ro/2, 0.0, 0.0] )
```
"""
function selectnode(fens::FENodeSet; kwargs...)
    nodelist = vselect(fens.xyz; kwargs...)
    nodelist = dropdims(reshape(nodelist,1,length(nodelist)), dims=1);
    return nodelist
end

"""
    selectelem(fens::FENodeSet, fes::T; kwargs...) where {T<:AbstractFESet}

Select finite elements.

# Arguments
- `fens` = finite element node set
- `fes` = finite element set
- `kwargs` = keyword arguments to specify the selection criteria

## Selection criteria

### facing
Select all "boundary" elements that "face" a certain direction.
```
exteriorbfl = selectelem(fens, bdryfes, facing=true, direction=[1.0, 1.0, 0.0]);
```
or
```
exteriorbfl = selectelem(fens, bdryfes, facing=true, direction=dout, dotmin = 0.99);
```
where
```
function dout(xyz)
    return xyz/norm(xyz)
end
```
and `xyz` is the location of the centroid  of  a boundary element.
Here the finite element is considered "facing" in the given direction if the dot
product of its normal and the direction vector is greater than `dotmin`.
The default value for `dotmin` is 0.01 (this corresponds to  almost 90 degrees
between the normal to the finite element  and the given direction).

This selection method makes sense only for elements that are  surface-like (i. e.
for boundary mmeshes).

### label
Select elements based on their label.
```
rl1 = selectelem(fens, fes, label=1)
```

### box, distance
Select elements based on some criteria that their nodes satisfy.  See the
function `selectnode()`.

Example:
Select all  elements whose nodes are closer than `R+inflate` from the point `from`.
```
linner = selectelem(fens, bfes, distance = R, from = [0.0 0.0 0.0],
  inflate = tolerance)
```

Example:

```
exteriorbfl = selectelem(fens, bdryfes,
   box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
```

### flood
Select all FEs connected together, starting from a given node.

Example:
Select all FEs connected together (Starting from node 13):
```julia
l = selectelem(fens, fes, flood = true, startnode = 13)
```

### Optional keyword arguments
Should we consider the element only if all its nodes are in?
- `allin` = Boolean: if true, then all nodes of an element must satisfy
the criterion; otherwise  one is enough.

# Output
`felist` = list of finite elements from the set that satisfy the criteria
"""
function selectelem(fens::FENodeSet, fes::T; kwargs...) where {T<:AbstractFESet}

    # smoothpatch
    #
    # Select all FEs that are part of a smooth surface.
    # For instance, starting from the finite element number 13, select all finite elements
    # whose normals defer from the normal of the Neighbor element by less than
    # 0.05 in the sense that dot(n1,n2)>sqrt(1-0.05^2).
    #     fe_select(fens,fes,struct ('flood', true, 'startfe', 13,  'normaldelta',0.05))
    # # Select all FEs "facing" in the direction [x(1),x(2),0] (away from the z-axis):
    #     fe_select(fens,fes,struct ('facing', true, 'direction', @(x)(+[x(1:2),0])))
    # Here x is the centroid of the nodes of each selected fe.
    # # Select all FEs "facing" in the direction [x(1),x(2),0] (away from the z-axis):
    #     fe_select(fens,fes,struct ('facing',1,'direction',@(x)(+[x(1:2),0]),'tolerance',0.01))
    # Here the fe is considered facing in the given direction if the dot
    # product of its normal and the direction vector is greater than tolerance.

    # Extract arguments
    allin = nothing; flood = nothing; facing = nothing; label = nothing;
    # nearestto = nothing; smoothpatch = nothing;
    startnode = 0; dotmin = 0.01
    overlappingbox = nothing
    for apair in pairs(kwargs)
        sy, val = apair
        if sy == :flood
            flood = val
        elseif sy == :facing
            facing = val
        elseif sy == :label
            label = val
        elseif sy == :overlappingbox
            overlappingbox = val
        elseif sy == :nearestto
            nearestto = val
        elseif sy == :allin
            allin = val
        end
    end

    if flood != nothing
        for apair in pairs(kwargs)
            sy, val = apair
            if sy == :startnode
                startnode = val
            end
        end
    end

    if facing != nothing
        facing = true;
        direction = nothing
        dotmin = 0.01;
        for apair in pairs(kwargs)
            sy, val = apair
            if sy == :direction
                direction = val
            elseif (sy == :dotmin) || (sy == :tolerance)# allow for obsolete keyword  to work
                dotmin = val
            end
        end
    end

    # if isfield(options,'smoothpatch')
    #     smoothpatch = true;
    #     if isfield(options,'normaldelta')
    #         normaldelta = options.normaldelta;
    #     else
    #         normaldelta= 0.05;
    #     end
    #     if isfield(options,'startfe')
    #         startfe = options.startfe;
    #     else
    #         error('Need the number of the starting finite element');
    #     end
    # end

    # if overlappingbox != nothing
    #     overlapping_box = true;
    #     # box=options.overlapping_box;
    #     # bounding_boxes =[];# precomputed bounding boxes of fes can be passed in
    #     # if isfield(options,'bounding_boxes')
    #     #     bounding_boxes = options.bounding_boxes;
    #     # end
    #     # inflate = 0;
    #     # if isfield(options,'inflate')
    #     #     inflate = (options.inflate);
    #     # end
    #     # box = inflate_box (box, inflate);
    # end

    # if isfield(options, 'nearestto')
    #     nearestto = true;
    #     locations = options.nearestto;
    # end

    # The  elements of this array are flipped from zero  when the element satisfies
    # the search condition.. This list is  eventually purged of the zero elements and
    # returned.
    felist = zeros(FInt,length(fes.conn));

    #     Select based on fe label
    if label!= nothing
        @assert length(fes.label) == length(fes.conn)
        for i=1:length(fes.conn)
            if label==fes.label[i]
                felist[i] =i;   # matched this element
            end
        end
        return  felist[findall(x->x!=0, felist)]; # return the nonzero element numbers
    end

    # Select by flooding
    if flood != nothing && (flood)
        @assert startnode > 0
        fen2fe = FENodeToFEMap(connasarray(fes), count(fens))
        felist = zeros(FInt, count(fes));
        pfelist = zeros(FInt, count(fes));
        felist[fen2fe.map[startnode]] .= 1;
        while true
            copyto!(pfelist, felist);
            markedl = findall(x -> x != 0, felist)
            for j = markedl
                for k = fes.conn[j]
                    felist[fen2fe.map[k]] .= 1;
                end
            end
            if sum(pfelist-felist) == 0 # If there are no more changes in this pass, we are done
                break;
            end
        end
        return findall(x -> x != 0, felist); # return the nonzero element numbers;
    end

    # Helper function: calculate the normal to a boundary finite element
    function normal(Tangents)
        sdim, mdim = size(Tangents);
        @assert (mdim==1) || (mdim==2)
        if     mdim==1 # 1-D fe
            N = [Tangents[2,1],-Tangents[1,1]];
            return N/norm(N);
        elseif mdim==2 # 2-D fe
            N= cross(Tangents[:,1],Tangents[:,2])
            return N/norm(N);
        end
    end

    # Select by in which direction the normal of the fes face
    if (facing != nothing) && (facing)
        xs =fens.xyz;
        sd =spacedim(fens);
        md =manifdim(fes);
        @assert (md == sd-1) "'Facing': only for Manifold dim. == Space dim.-1"
        param_coords =zeros(FFlt,1,md);
        Need_Evaluation = (typeof(direction) <: Function);
        if ( !Need_Evaluation)
            d = reshape(direction,1,sd)/norm(direction);
        end
        Nder = bfundpar(fes, vec(param_coords));
        for i=1:length(fes.conn)
            xyz =xs[[j for j in fes.conn[i]],:];
            Tangents =xyz'*Nder;
            N = normal(Tangents);
            if (Need_Evaluation)
                d = direction(mean(xyz, dims = 1));
                d = reshape(d,1,sd)/norm(d);
            end
            if (dot(vec(N),vec(d)) > dotmin)
                felist[i]=i;
            end
        end
        return  felist[findall(x->x!=0, felist)]; # return the nonzero element numbers
    end


    #   # Select by the change in normal
    #   if (smoothpatch)
    #     mincos=sqrt(1-normaldelta^2);
    #     xs=fens.xyz;
    #     sdim =size(xs,2);
    #     mdim=fes.dim;
    #     if (mdim~=sdim-1)
    #       error ('"Smoothpatch" selection of fes make sense only for Manifold dimension ==Space dimension-1')
    #     end
    #     fen2fe_map=fenode_to_fe_map (struct ('fes',fes));
    #     gmap=fen2fe_map.map;
    #     conn=fes.conn;
    #     param_coords =zeros(1,mdim);# This is a hack: this may not be the proper location for all element types
    #     Nder = bfundpar (fes, param_coords);
    #     felist= zeros(1, size(conn,1));
    #     felist(startfe)=1;
    #     while true
    #       pfelist=felist;
    #       markedl=find(felist~=0);
    #       for j=markedl
    #         xyz =xs(conn(j,:),:);
    #         Tangents =xyz'*Nder;
    #         Nj = normal (Tangents,sdim, mdim);
    #         for k=conn(j,:)
    #           for ke=gmap{k}
    #             xyz =xs(conn(ke,:),:);
    #             Tangents =xyz'*Nder;
    #             Nk = normal (Tangents,sdim, mdim);
    #             if (dot(Nj,Nk)>mincos)
    #               felist(ke)=1;
    #             end
    #           end
    #         end
    #       end
    #       if sum(pfelist-felist)==0, break; end
    #     end
    #     felist =find(felist~=0);
    #     return;
    #   end


    # Select all FEs whose bounding box overlaps given box
    if (overlappingbox != nothing)
        bbox = zeros(2 * size(fens.xyz, 2))
        for i = 1:length(fes.conn)
            bbox = initbox!(bbox, vec(fens.xyz[fes.conn[i][1], :]))
            for k = 2:length(fes.conn[i])
                bbox = updatebox!(bbox, fens.xyz[fes.conn[i][k], :])
            end
            if boxesoverlap(overlappingbox, bbox)
                felist[i]=i; # this element overlaps the box
            end
        end
        return  felist[findall(x->x!=0, felist)]; # return the nonzero element numbers
    end


    #   # get the fe nearest to the supplied point
    #   if (nearestto)

    #     felist = [];

    #     # get the locations of all the nodes
    #     xs = fens.xyz;
    #     #     Disqualify all the nodes that are not connected to the finite
    #     #     elements on input
    #     cn = connected_nodes(fes);
    #     xs(setdiff(1:size(xs,1),cn),:) =Inf;

    #     # get the connectivity of all the FEs
    #     conns = fes.conn;

    #     for i = 1:size(locations,1)

    #       # Get the smallest distances between the node locations and the
    #       # desired location using ipdm (by John D'Errico)
    #       distance = ipdm(xs, locations(i,:));

    #       # Get array index of smallest distance to location from
    #       # all nodes
    #       [junk,IX] = min(distance);

    #       # Find the fes connected to the nearest node
    #       [crows, ccols] = find(conns == IX(1));

    #       if numel(crows) == 1
    #         # if there's only one, this is superb, we can add it to the
    #         # list and continue to the next location
    #         felist(end+1) = crows;
    #       else

    #         # otherwise we must determine which cell is closest based on
    #         # the other connected nodes
    #         nearconns = conns(crows, :);

    #         nearconndist = zeros(size(nearconns));

    #         # Get the distances between all the nodes in the nearest
    #         # FEs and the nearest node to the location
    #         for j = 1:size(nearconns, 1)
    #           for k = 1:size(nearconns, 2)
    #             nearconndist(j,k) = ipdm(xs(nearconns(j,k), :), locations(i,:));
    #           end
    #         end

    #         # set the distance from the nearest node to in the cell to
    #         # the location to infinity
    #         nearconndist(nearconns == IX(1)) = inf;

    #         # order the FEs by the smallest distance of any connected
    #         # node from the nearest node
    #         [junk,IX] = sort(min(nearconndist,[],2));

    #         # return the index in conns to the first of these ordered
    #         # FEs
    #         if (~isempty(crows))
    #           felist(end+1) = crows(IX(1));
    #         end

    #       end

    #     end

    #     return;

    #   end

    #  Default:   Select based on location of nodes
    #   Should we consider the element only if all its nodes are in?
    allinvalue = (allin == nothing) || ((allin != nothing) && (allin))
    # Select elements whose nodes are in the selected node list
    nodelist = selectnode(fens; kwargs...);
    nper = nodesperelem(fes)
    for i = 1:length(fes.conn)
        found = 0
        for nd in fes.conn[i]
            if nd in nodelist
                found = found + 1
            end
        end
        if allinvalue
            if found == nper
                felist[i] =i;
            end
        else
            if found >= 1
                felist[i] =i;
            end
        end
    end
    return felist[findall(x->x!=0, felist)]; # return the nonzero element numbers

end

function _compute_vlist!(vlist::Vector{FInt}, abox::T, sdim::FInt, v::FFltMat) where {T}
    # Helper functions
    @inline inrange(rangelo,rangehi,x) = (rangelo <= x <= rangehi)
    nn = 0
    for j in 1:size(v,1)
        matches = true
        for i in 1:sdim
            if !inrange(abox[2*i-1], abox[2*i], v[j, i])
                matches = false; break
            end
        end
        if matches
            nn = nn + 1; vlist[nn] = j;
        end
    end
    return vlist, nn
end

"""
    vselect(v::FFltMat; kwargs...)

Select locations (vertices) from the array based on some criterion.

See the function `selectnode()` for examples of the criteria that can be
used to search vertices.
"""
function vselect(v::FFltMat; kwargs...)




    # Extract arguments
    box = nothing; distance = nothing; from = nothing; plane  =  nothing;
    thickness = nothing; nearestto = nothing; inflate = 0.0;
    for apair in pairs(kwargs)
        sy, val = apair
        if sy == :box
            box = val
        elseif sy == :distance
            distance = val
        elseif sy == :from
            from = val
        elseif sy == :plane
            plane = val
        elseif sy == :thickness
            thickness = val
        elseif sy == :nearestto
            nearestto = val
        elseif sy == :inflate
            inflate = val
        end
    end

    # Did we get an inflate value
    inflatevalue = 0.0;
    if inflate != nothing
        inflatevalue = FFlt(inflate);
    end

    # Initialize the output list
    vlist = zeros(FInt, size(v,1)); nn = 0;


    # Process the different options
    if box != nothing
        sdim = size(v,2)
        dim = Int(round(length(box)/2.))::FInt;
        @assert dim == sdim "Dimension of box not matched to dimension of array of vertices"
        abox = vec(box)::FFltVec
        inflatebox!(abox, inflatevalue)
        vlist, nn = _compute_vlist!(vlist, abox, sdim, v)
    elseif distance != nothing
        fromvalue =fill!(deepcopy(v[1,:]), 0.0);
        if from!=nothing
            fromvalue = from;
        end
        fromvalue = reshape(fromvalue, (size(v[1,:])))
        d=distance+inflatevalue;
        for i=1:size(v,1)
            if norm(fromvalue-v[i,:])<d
                nn =nn +1; vlist[nn] =i;
            end
        end
    elseif plane != nothing
        normal = plane[1:end-1];
        normal = vec(normal/norm(normal));
        thicknessvalue = 0.0;
        if thickness != nothing
            thicknessvalue = thickness;
        end
        t = thicknessvalue+inflatevalue;
        distance = plane[end];
        for i = 1:size(v,1)
            ad = dot(v[i,:], normal);
            if abs(distance-ad)<t
                nn = nn +1; vlist[nn] =i;
            end
        end
    elseif nearestto != nothing
        location = vec(nearestto);
        distance =  zeros(size(v,1));
        for i=1:size(v,1)
            distance[i] = norm(location-vec(v[i,:]))
        end
        Mv,j = findmin(distance)
        vlist[1] = j;
        nn=1;
    end
    if (nn==0)
        vlist = FInt[];# nothing matched
    else
        vlist = vlist[1:nn];
    end
    return vlist
end

"""
    findunconnnodes(fens::FENodeSet, fes::AbstractFESet)

Find nodes that are not connected to any finite element.

connected = array is returned which is for the node k either true (node k is
     connected), or false (node k is not connected).

Let us say there are nodes not connected to any finite element that you
would like to remove from the mesh: here is how that would be
accomplished.
"""
function findunconnnodes(fens::FENodeSet, fes::AbstractFESet)
    connected = trues(count(fens));
    fen2fem = FENodeToFEMap(connasarray(fes), count(fens))
    for i=1:length(fen2fem.map),
        connected[i] = (!isempty(fen2fem.map[i]));
    end
    return connected
end


end

# NOTE: This operation is probably best done with the  node-to-element  map.
# """
#     connectedelems(fes::AbstractFESet, node_list::FIntVec)
#
# Extract the numbers of the finite elements connected to given nodes.
#
# Extract the list of numbers for the fes  that are connected to given nodes.
#
# Warning: this tends to be an expensive operation.
# """
# function connectedelems(fes::AbstractFESet, node_list::FIntVec)
#   cg=zeros(FInt,size(fes.conn,1));
#   for j=1:size(fes.conn,1)
#     cg[j]= length( intersect(fes.conn[j,:], node_list) );
#   end
#   cg =find_n(cg);
# end
