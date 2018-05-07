"""
Module for remeshing tetrahedral triangulations.
"""
module TetRemeshingModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import Base.length
import Base.push!
import Base.getindex
import Base.copy!
import FinEtools.FENodeToFEMapModule: FENodeToFEMap
import FinEtools.MeshTetrahedronModule: T4meshedges
import FinEtools.MeshModificationModule: interior2boundary
import LinearAlgebra: norm

"""
    _IntegerBuffer

This data structure is introduced in order to avoid constant allocation/deallocation of integer  arrays.  
The buffer is allocated once and then reused  millions of times. The memory is never  released.
The buffer needs to be accessed  through the methods  below, not directly.
"""
mutable struct _IntegerBuffer
    a::Array{Int,1}
    pa::Int
end

function push!(b::_IntegerBuffer, i::Int)
    b.pa = b.pa + 1 
    b.a[b.pa] = i 
    return b
end

function empty!(b::_IntegerBuffer)
    b.pa = 0
    return b
end

function trim!(b::_IntegerBuffer, i::Int)
    b.pa = i
    return b
end

length(b::_IntegerBuffer) = b.pa

function getindex(b::_IntegerBuffer, i::Int)
    @assert i <= length(b)
    return b.a[i]
end 

function copyto!(d::_IntegerBuffer, s::_IntegerBuffer)
    @assert length(d.a) == length(s.a)
    empty!(d)
    @inbounds for i = 1:length(s) # @inbounds 
        push!(d, s.a[i])
    end
    return d
end

function tetvtimes6(vi1::FFltVec, vi2::FFltVec, vi3::FFltVec, vi4::FFltVec)
    # local one6th = 1.0/6
    # @assert size(X, 1) == 4
    # @assert size(X, 2) == 3
    @inbounds let
        A1 = vi2[1]-vi1[1]; 
        A2 = vi2[2]-vi1[2]; 
        A3 = vi2[3]-vi1[3]; 
        B1 = vi3[1]-vi1[1]; 
        B2 = vi3[2]-vi1[2]; 
        B3 = vi3[3]-vi1[3]; 
        C1 = vi4[1]-vi1[1]; 
        C2 = vi4[2]-vi1[2]; 
        C3 = vi4[3]-vi1[3]; 
        return ((-A3*B2+A2*B3)*C1 +  (A3*B1-A1*B3)*C2 + (-A2*B1+A1*B2)*C3);
    end
end

function checkvolumes(when, vt, t)
    for i = 1:size(t, 1)
        if any(x->x==0, t[i,:])
        else 
            if tetvtimes6(vt[:,t[i,1]], vt[:,t[i,2]], vt[:,t[i,3]], vt[:,t[i,4]]) < 0.0
                error("*** $when Negative volume = $([t[i,:]])")
            end
        end 
    end
end

"""
    coarsen(t::Array{Int, 2}, v::Array{Float64, 2}, tmid::Vector{Int}, options)

Coarsen a T4 (tetrahedral) mesh.  

## Arguments
t = array with vertex numbers, one per tetrahedron in each row
v= array of vertex coordinates, one per vertex in each row
tmid = tetrahedron material identifiers, one for tetrahedron in each row
options = structure with optional fields
bv=array for all vertices in the input array v.  
       Is the vertex on the boundary? True or false. The boundary vertices 
       are in layer 1. The vertex layer numbers increase away from 
       the boundary,  and the bigger the vertex layer number the bigger 
       the mesh size.
desired_ts=mesh size desired for vertices in layer 1, here mesh size is 
       the length of an edge
stretch=the mesh size should increase by this much within one layer of  
       elements, default is 1.25
nblayer=number of boundary layers where the mesh size should not increase, 
       default is 1
surface_coarsening = is the coarsening intended for the interior or for 
       the surface? default is false, which means the default 
       is coarsening of the interior.
preserve_thin= when coarsening the surface, should features which are  thin 
       be unaffected by the coarsening? Here thin means consisting of 
       only "surface" vertices.
vertex_weight= weight of vertices, one per  vertex; weight <= 1.0 is ignored, 
       but weight>1.0  is used to "stretch" edges connected  to the vertex.
       In this way one may preserve certain edges by assigning  larger 
       weight to their vertices. Default is vertex_weight= [] (which means 
       ignore the vertex weights)

## Output
t,v,tmid = new arrays
"""
function coarsen(t::Array{Int, 2}, inputv::Array{Float64, 2}, tmid::Vector{Int}; bv::Vector{Bool} = Bool[], desired_ts::Number = 0.0, stretch::Number = 1.25, nblayer::Int = 1, surface_coarsening::Bool = false, preserve_thin::Bool = false, vertex_weight::Vector{Float64} = Float64[], reportprogress::F = n -> nothing) where {F<:Function}
    vt = deepcopy(transpose(inputv))  # Better locality of data can be achieved with vertex coordinates in columns
    nv = size(inputv, 1)
    if length(bv) == nv
        vlayer = Int[i == true ? 1 : 0 for i in bv]
    else 
        vlayer = Int[]
    end 
    if length(vertex_weight) != nv
        vertex_weight = ones(nv)
    end 
    m = FENodeToFEMap(t, nv)
    v2t = deepcopy(m.map);# Map vertex to tetrahedron
    e = T4meshedges(t);
    m = FENodeToFEMap(e, nv)
    v2e = deepcopy(m.map);# Map vertex to edge
    
    # Figure out the  vertex layer numbers which are the guides to coarsening
    if (isempty(vlayer))
        # # Extract the boundary faces
        f = interior2boundary(t, [1 3 2; 1 2 4; 2 3 4; 1 4 3])
        vlayer = zeros(Int, nv);
        if (surface_coarsening)
            i = setdiff(collect(1:nv),vec(f));
            vlayer[i] .= 1;# Mark all vertices in the interior (layer 1) 
            vlayer[f] .= 2;# Mark all vertices on the boundary (layer 2) 
            if (preserve_thin) # thin features should be preserved, 
                # mark them as interior (layer 1) so the algorithm does not touch them
                oldvlayer = deepcopy(vlayer);
                for j = vec(f)
                    # is any vertex connected to j  in the interior?
                    connected =false;
                    for k=1:length(v2e[j])
                        if (oldvlayer(e[v2e[j][k],1]) == 1) 
                            connected = true; break;
                        end
                        if (oldvlayer(e[v2e[j][k],2]) == 1) 
                            connected = true; break;
                        end
                    end
                    if (!connected)
                        vlayer[j]=1; # mark it as interior to prevent coarsening
                    end
                end
            end
        else # not (surface_coarsening)
            vlayer[f] .= 1;# Mark all vertices on the boundary (layer 1) 
            # The remaining (interior) vertices will be numbered below 
            # by  distance from the boundary 
        end
    end
    
    # Compute the vertex layer numbers for vertices for which they have not 
    # been supplied (which is marked by 0).
    currvlayer = 1;
    while true
        # mark layer next to current boundary
        oldvlayer = deepcopy(vlayer);
        for i=1:size(vlayer,1)
            if oldvlayer[i] == currvlayer
                for j=1:length(v2t[i])
                    k=v2t[i][j];
                    for m=1:4
                        if (oldvlayer[t[k,m]]==0)
                            vlayer[t[k,m]] = currvlayer+1;
                        end
                    end
                end
            end
        end
        if isempty(findall(p->p==0, vlayer)) || (norm(vlayer-oldvlayer)==0)
            break;
        end
        currvlayer =currvlayer +1;
    end
    currvlayer = currvlayer +1;
    # Compute the desired mesh size at a given layer
    desiredes = zeros(currvlayer);
    desiredes[1] = desired_ts;
    for layer =2:currvlayer
        s=0; finalr = 0
        for r=1:currvlayer
            s=s+stretch^r;
            finalr = r
            if (s>= layer-1)
                break
            end
        end
        finalr=finalr-1;
        desiredes[layer] = desired_ts*stretch^finalr;
    end
    
    # Initialize edge lengths, edge vertex counts
    es = zeros(size(e,1));
    elayer = zeros(Int, size(e,1));
    es = edgelengths(collect(1:size(e,1)), vt, e, vertex_weight);
    elayer = vec(minimum(hcat(vlayer[e[:,1]],vlayer[e[:,2]]), dims = 2))
    everfailed = zeros(Bool, size(e,1));
    minne = 200;    maxne = 400;# 36
    availe = _IntegerBuffer(zeros(Int, size(e,1)), 0) # Flexible  buffer to avoid allocations
    selist = _IntegerBuffer(zeros(Int, size(e,1)), 0) # Flexible  buffer to avoid allocations
    
    checkvolumes("Before reduction starts", vt, t)

    # Let's get down to business
    previouscurrvlayer = 0
    pass =1;
    while true
        availe, currvlayer = availelist!(availe, selist, elayer, currvlayer, minne, maxne, es, nblayer, desiredes, everfailed);
        if currvlayer != previouscurrvlayer
            reportprogress(currvlayer)
        end 
        if (length(availe)==0)  
            break; # Done. Hallelujah!
        end 
        while true
            selist = sortedelist!(selist, elayer, es, desiredes, availe);
            Change =  false; 
            for i=1:length(selist)
                if (collapseedge!(e, es, elayer, vlayer, t, vt, v2t, v2e, vertex_weight, selist[i]))
                    Change = true;  break; 
                end
                everfailed[selist[i]] = true; # the collapse failed for this edge
            end
            if  (!Change)
                break;
            end
        end
        pass = pass + 1; previouscurrvlayer = currvlayer
    end
    reportprogress(0) # we are done
    # Note  that we are reverting the transpose of the vertex array here
    checkvolumes("Before clean", vt, t)
    return cleanoutput(t,deepcopy(transpose(vt)),tmid);
end

function collapseedge!(e::Array{Int, 2}, es::Vector{Float64}, elayer::Vector{Int}, vlayer::Vector{Int}, t::Array{Int, 2}, vt::AbstractArray{Float64, 2}, v2t::Vector{Vector{FInt}}, v2e::Vector{Vector{FInt}}, vertex_weight::Vector{Float64}, dei::Int)
    result = false;
    # if the operation would result in inverted tetrahedra, cancel it
    de1, de2 = e[dei, 1], e[dei, 2];
    if anynegvol1(t, v2t[de2], vt, de2, de1)
        de1, de2 = e[dei, 2], e[dei, 1];;# Try the edge the other way
        if anynegvol1(t, v2t[de2], vt, de2, de1)
            return false; # the collapse failed
        end
    end
    vi1 = de1; vi2 = de2; # the kept Vertex,  and the replaced Vertex
    # Modify t: switch the references to the replaced vertex vi2
    mtl = unique(vcat(v2t[vi1],v2t[vi2]));
    for k = 1:length(mtl)
        for i = 1:4
            if t[mtl[k], i] == vi2
                t[mtl[k], i] = vi1; 
            end 
        end
    end
    # Modify e: switch the references to the replaced vertex vi2
    mel = unique(vcat(v2e[vi1],v2e[vi2]));
    for k=1:length(mel)
        for i = 1:2
            if e[mel[k], i] == vi2
                e[mel[k], i] = vi1; 
            end 
        end
    end
    # Delete the collapsed tetrahedra from the vertex-to-tet map
    dtl = intersect(v2t[vi1],v2t[vi2]); # list of tets that connect the verts on the collapsed edge
    vl = unique(vec(t[dtl,:])); # vertices incident on the collapsed tetrahedra
    for i=1:length(vl) # Delete the collapsed tetrahedra
        filter!(p -> !(p in dtl), v2t[vl[i]])
    end
    # Delete edges which are merged by the collapse
    del = v2e[vi2]; # vi2 is the vertex that is to be deleted
    for k =1:length(del)
        i = del[k];
        (e[i,1] == vi2) ? ov = e[i,2] : ov = e[i,1]; 
        if ov in vl
            e[i,:] .= 0;# mark as deleted
            es[i] = Inf;# indicate the deleted edge
            elayer[i] = 0;# indicate the deleted edge
        end
    end
    t[dtl,:] .= 0;# Mark deleted tetrahedra
    e[dei,:] .= 0;# Mark deleted edge
    # Update the vertex-2-tet  map
    v2t[vi1]= setdiff(mtl, dtl);
    v2t[vi2]= Int[];# this vertex is gone
    # Update the vertex-2-edge  map
    v2e[vi1]= setdiff(mel, dei);
    v2e[vi2]= Int[];# this vertex is gone
    #     v(vi1,:) = nv1;# new vertex location
    vt[:,vi2] .= Inf;# Indicate invalid vertex
    # update edge lengths
    for k=1:length(v2e[vi1])
        i=v2e[vi1][k];
        if (e[i,1]==0)
            es[i] = Inf;# indicate the deleted edge
            elayer[i] = 0;
        else
            es[i] = elength(vt, vertex_weight, e[i,1], e[i,2])
            elayer[i] = minimum(vlayer[e[i,:]]);
        end
    end
    es[dei] = Inf;# indicate the deleted edge
    elayer[dei] = 0; # indicate the deleted edge
    return  true;
end

function edgelengths(ens::Vector{Int}, vt::AbstractArray{Float64, 2}, e::Array{Int, 2}, vertex_weight::Vector{Float64})
    eLengths = zeros(length(ens))
    for i = 1:length(ens)
        en = ens[i]
        eLengths[i] = elength(vt, vertex_weight, e[en,1], e[en,2])
    end
    return eLengths
end

# Weighted length of the edge between 2 vertices
function elength(vt::AbstractArray{Float64, 2}, vertex_weight::Vector{Float64}, i1::Int, i2::Int)
    @inbounds return max(vertex_weight[i1], vertex_weight[i2]) * 
        sqrt((vt[1,i2] - vt[1,i1])^2 + (vt[2,i2] - vt[2,i1])^2 + (vt[3,i2] - vt[3,i1])^2);
end

function sortedelist!(selist::_IntegerBuffer, elayer::Vector{Int}, es::Vector{Float64}, desiredes::Vector{Float64}, availe::_IntegerBuffer)
    empty!(selist)
    for i11=1:length(availe)
        k11=availe[i11];
        if (elayer[k11] >0)
            if  es[k11] < desiredes[elayer[k11]]
                push!(selist, k11);
            end
        end
    end
    return selist;
end

function availelist!(availe::_IntegerBuffer, selist::_IntegerBuffer, elayer::Vector{Int}, currvlayer::Int, minnt::Int, maxnt::Int, es::Vector{Float64}, nblayer::Int, desiredes::Vector{Float64}, everfailed::Array{Bool,1})
    newcurrvlayer = 0
    empty!(selist)
    for layer = currvlayer:-1:nblayer+1 # This can be changed to allow for more or less coarsening
        empty!(availe)
        @inbounds for i = 1:length(elayer) # @inbounds 
            if (layer <= elayer[i]) && (!everfailed[i])
                push!(availe, i)
            end 
        end
        selist = sortedelist!(selist, elayer, es, desiredes, availe);
        newcurrvlayer = layer;
        if (length(selist) >= minnt)
            break;
        end
    end
    return copyto!(availe, trim!(selist, min(length(selist), maxnt))), newcurrvlayer
end

function anynegvol1(t::Array{Int,2}, whichtets::Vector{Int}, vt::AbstractArray{Float64,2}, whichv::Int, otherv::Int)
    for iS1 in whichtets
        i1, i2, i3, i4 = t[iS1,:]; # nodes of the tetrahedron
        if (i1 == whichv) 
            i1 = otherv
        end
        if (i2 == whichv) 
            i2 = otherv
        end
        if (i3 == whichv) 
            i3 = otherv
        end
        if (i4 == whichv) 
            i4 = otherv
        end
        if tetvtimes6(vt[:,i1], vt[:,i2], vt[:,i3], vt[:,i4]) < 0.0
            return true
        end 
    end
    return false;
end

function cleanoutput(t::Array{Int,2}, v::Array{Float64,2}, tmid::Array{Int,1})
    nn = zeros(Int, size(v,1));
    nv = deepcopy(v);
    k=1;
    for i=1:size(v,1)
        if (v[i,1] != Inf) # Is this an active vertex?
            nn[i] = k;
            nv[k,:] =v[i,:];
            k=k+1;
        end
    end
    nnv=k-1;
    v = nv[1:nnv,:];
    # delete connectivities of collapsed tetrahedra, and renumber nodes
    nt = deepcopy(t); fill!(nt, 0)
    ntmid = deepcopy(tmid); fill!(ntmid, 0)
    k=1;
    for i=1:size(t,1)
        if (t[i,1] != 0)# not deleted
            if (!isempty(findall(p -> p == 0, nn[t[i,:]])))
                # error('Referring to deleted vertex')
                t[i,:] = 0;
            else
                j =nn[t[i,:]]
                nt[k,:]=j;
                ntmid[k]=tmid[i];
                k=k+1;
            end
        end
    end
    t = nt[1:k-1,:];
    # delete unconnected vertices
    uv = unique(vec(t));
    if (length(uv) != size(v,1)) # there may be unconnected vertices
        nn = zeros(Int, size(v,1),1);
        nn[uv] = collect(1:length(uv));
        nv = deepcopy(v);
        for i=1:size(v,1)
            if (nn[i] != 0)
                nv[nn[i],:] =v[i,:];
            end
        end
        v = nv[1:length(uv),:];
        for i=1:size(t,1)
            t[i,:] = nn[t[i,:]];
        end
    end
    tmid = ntmid[1:k-1];
    return t, v, tmid
end

end # module

# function t4_e2(t::Array{Int, 2})
#     ec = [  1  2
#             2  3
#             3  1
#             4  1
#             4  2
#             4  3];
#     e = vcat(t[:,ec[1,:]], t[:,ec[2,:]], t[:,ec[3,:]], t[:,ec[4,:]], t[:,ec[5,:]], t[:,ec[6,:]])
#     e = sort(e, 2);
#     ix = sortperm(e[:,1]);
#     e = e[ix,:];
#     ue = deepcopy(e)
#     i = 1;
#     n=1;
#     while n <= size(e,1)
#         c = ue[n,1];
#         m = n+1;
#         while m <= size(e,1)
#             if (ue[m,1] != c)
#                 break; 
#             end
#             m = m+1;
#         end
#         us = unique(ue[n:m-1,2], 1);# the transcription below is a lot faster
#         ls =length(us);
#         e[i:i+ls-1,1] = c;
#         e[i:i+ls-1,2] = s;
#         i = i+ls;
#         n = m;
#     end
#     e = e[1:i-1,:];
# end

# """
# anynegvol(t, v, whichv, v1)

# Check that the new location `v1` of the vertex `whichv` does not result 
# in a negative volume  for any of the tetrahedra connected to `whichv`.  

# This is a heavily used function, and hence speed and absence of
# memory allocation is at a premium.
# """
# function anynegvol(t::Array{Int,2}, whichtets::Vector{Int}, vt::Array{Float64,2}, whichv::Int, otherv::Int)
# i1, i2, i3, i4 = 0, 0, 0, 0
# Volume6 = 0.0 # @inbounds Volume6 = let
# for iS1 in whichtets
#     i1, i2, i3, i4 = t[iS1,:]; # nodes of the tetrahedron
#     if (i1 == whichv) 
#         @inbounds Volume6 = let
#             A1 = vt[1,i2]-vt[1,otherv]; 
#             A2 = vt[2,i2]-vt[2,otherv]; 
#             A3 = vt[3,i2]-vt[3,otherv]; 
#             B1 = vt[1,i3]-vt[1,otherv]; 
#             B2 = vt[2,i3]-vt[2,otherv]; 
#             B3 = vt[3,i3]-vt[3,otherv]; 
#             C1 = vt[1,i4]-vt[1,otherv]; 
#             C2 = vt[2,i4]-vt[2,otherv]; 
#             C3 = vt[3,i4]-vt[3,otherv]; 
#             ((-A3*B2+A2*B3)*C1 +  (A3*B1-A1*B3)*C2 + (-A2*B1+A1*B2)*C3)
#         end
#     end
#     if (i2 == whichv) 
#         @inbounds Volume6 = let
#             A1 = vt[1,otherv]-vt[1,i1]; 
#             A2 = vt[2,otherv]-vt[2,i1]; 
#             A3 = vt[3,otherv]-vt[3,i1]; 
#             B1 = vt[1,i3]-vt[1,i1]; 
#             B2 = vt[2,i3]-vt[2,i1]; 
#             B3 = vt[3,i3]-vt[3,i1]; 
#             C1 = vt[1,i4]-vt[1,i1]; 
#             C2 = vt[2,i4]-vt[2,i1]; 
#             C3 = vt[3,i4]-vt[3,i1]; 
#             ((-A3*B2+A2*B3)*C1 +  (A3*B1-A1*B3)*C2 + (-A2*B1+A1*B2)*C3)
#         end
#     end
#     if (i3 == whichv) 
#         @inbounds Volume6 = let
#             A1 = vt[1,i2]-vt[1,i1]; 
#             A2 = vt[2,i2]-vt[2,i1]; 
#             A3 = vt[3,i2]-vt[3,i1]; 
#             B1 = vt[1,otherv]-vt[1,i1]; 
#             B2 = vt[2,otherv]-vt[2,i1]; 
#             B3 = vt[3,otherv]-vt[3,i1]; 
#             C1 = vt[1,i4]-vt[1,i1]; 
#             C2 = vt[2,i4]-vt[2,i1]; 
#             C3 = vt[3,i4]-vt[3,i1]; 
#             ((-A3*B2+A2*B3)*C1 +  (A3*B1-A1*B3)*C2 + (-A2*B1+A1*B2)*C3)
#         end
#     end
#     if (i4 == whichv) 
#         @inbounds Volume6 = let
#             A1 = vt[1,i2]-vt[1,i1]; 
#             A2 = vt[2,i2]-vt[2,i1]; 
#             A3 = vt[3,i2]-vt[3,i1]; 
#             B1 = vt[1,i3]-vt[1,i1]; 
#             B2 = vt[2,i3]-vt[2,i1]; 
#             B3 = vt[3,i3]-vt[3,i1]; 
#             C1 = vt[1,otherv]-vt[1,i1]; 
#             C2 = vt[2,otherv]-vt[2,i1]; 
#             C3 = vt[3,otherv]-vt[3,i1]; 
#             ((-A3*B2+A2*B3)*C1 +  (A3*B1-A1*B3)*C2 + (-A2*B1+A1*B2)*C3)
#         end
#     end
#     if (Volume6 < 0)
#          return true;
#     end
# end
# return false;
# end


# function dumparray(a::Array{T,1}, afile) where {T}
#     fid=open(afile * ".csv","w");
#     if (fid==-1)
#         error("Could not open " * afile)
#         return nothing
#     end
#     for i= 1:size(a, 1)
#         print(fid,a[i],"\n");
#     end
#     fid=close(fid);
# end

# function dumparray(a::Array{T,2}, afile) where {T}
#     fid=open(afile * ".csv","w");
#     if (fid==-1)
#         error("Could not open " * afile)
#         return nothing
#     end
#     for i= 1:size(a, 1)
#         for j= 1:size(a,2)-1
#             print(fid,a[i,j],",");
#         end
#         print(fid,a[i,end],",\n");
#     end
#     fid=close(fid);
# end

# julia> include("src\\TetRemeshingModule.jl"); include("test/playground.jl")
# WARNING: replacing module TetRemeshingModule.
# Mesh size: initial = 3000000
# 26.25.24.23.22.21.20.19.18.17.16.15.14.13.12.11.10.9.8.7.6.5.4.
# 3.2.
# 187.307604 seconds (520.54 M allocations: 141.169 GiB, 41.79% gc time)
# Mesh size: final = 401937 [187.5880000591278 sec]

# WARNING: replacing module TetRemeshingModule.
# Mesh size: initial = 3000000
# 26.25.24.23.22.21.20.19.18.17.16.15.14.13.12.11.10.9.8.7.6.5.4.
# 3.2.
# Mesh size: final = 401937 [133.66499996185303 sec]
# V = 0.013499999999977705 compared to 0.0135
# Task (runnable) @0x000000000af5dfb0

# With optimization  3  (passing the whole t to anyneg)
# PetrKrysl@Firebolt MINGW64 ~/.julia/v0.7/FinEtools (master)
# $ /c/Users/PetrKrysl/AppData/Local/Julia-0.7.0-DEV/bin/julia.exe
#                _
#    _       _ _(_)_     |  A fresh approach to technical computing
#   (_)     | (_) (_)    |  Documentation: https://docs.julialang.org
#    _ _   _| |_  __ _   |  Type "?help" for help.
#   | | | | | | |/ _` |  |
#   | | |_| | | | (_| |  |  Version 0.7.0-DEV.2098 (2017-10-10 11:37 UTC)
#  _/ |\__'_|_|_|\__'_|  |  Commit 546a801260* (16 days old master)
# |__/                   |  x86_64-w64-mingw32

# julia> include("src\\TetRemeshingModule.jl"); include("test/playground.jl")
# Mesh size: initial = 3000000
# 26.25.24.23.22.21.20.19.18.17.16.15.14.13.12.11.10.9.8.7.6.5.4.
# 3.2.
# Mesh size: final = 401937 [107.39300012588501 sec]
# V = 0.013499999999977705 compared to 0.0135
# Task (runnable) @0x000000000bf4cb90

# julia> include("src\\TetRemeshingModule.jl"); include("test/playground.jl")
# WARNING: replacing module TetRemeshingModule.
# Mesh size: initial = 3000000
# 26.25.24.23.22.21.20.19.18.17.16.15.14.13.12.11.10.9.8.7.6.5.4.
# 3.2.
# Mesh size: final = 401937 [93.27800011634827 sec]
# V = 0.013499999999977705 compared to 0.0135
# Task (runnable) @0x000000000da705d0
