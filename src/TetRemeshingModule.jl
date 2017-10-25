module TetRemeshingModule

using FinEtools
using FinEtools.FENodeToFEMapModule
using FinEtools.MeshTetrahedronModule

function dumparray(a::Array{T,1}, afile) where {T}
    fid=open(afile * ".csv","w");
    if (fid==-1)
        error("Could not open " * afile)
        return nothing
    end
    for i= 1:size(a, 1)
        print(fid,a[i],"\n");
    end
    fid=close(fid);
end

function dumparray(a::Array{T,2}, afile) where {T}
    fid=open(afile * ".csv","w");
    if (fid==-1)
        error("Could not open " * afile)
        return nothing
    end
    for i= 1:size(a, 1)
        for j= 1:size(a,2)-1
            print(fid,a[i,j],",");
        end
        print(fid,a[i,end],",\n");
    end
    fid=close(fid);
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
function coarsen(t::Array{Int, 2}, v::Array{Float64, 2}, tmid::Vector{Int}; bv::Vector{Bool} = Bool[], desired_ts::Number = 0.0, stretch::Number = 1.25, nblayer::Int = 1, surface_coarsening::Bool = false, preserve_thin::Bool = false, vertex_weight::Vector{Float64} = Float64[])
    if length(bv) == size(v, 1)
        vlayer = [i == true ? 1 : 0 for i in bv]
    else 
        vlayer = []
    end 
    if length(vertex_weight) != size(v, 1)
        vertex_weight = ones(size(v, 1))
    end 
    m = FENodeToFEMap(t, size(v,1))
    v2t = deepcopy(m.map);# Map vertex to tetrahedron
    e = MeshTetrahedronModule.T4meshedges(t);
    m = FENodeToFEMap(e, size(v,1))
    v2e = deepcopy(m.map);# Map vertex to edge
    
    # Figure out the  vertex layer numbers which are the guides to coarsening
    if (isempty(vlayer))
        # # Extract the boundary faces
        f = interior2boundary(t, [1 3 2; 1 2 4; 2 3 4; 1 4 3])
        vlayer = zeros(Int, size(v,1));
        if (surface_coarsening)
            i = setdiff(collect(1:size(v,1)),vec(f));
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
            vlayer[f] = 1;# Mark all vertices on the boundary (layer 1) 
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
        if isempty(find(p->p==0, vlayer)) || (norm(vlayer-oldvlayer)==0)
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
        desiredes[layer] .= desired_ts*stretch^finalr;
    end
    
    # Initialize edge lengths, edge vertex counts
    es = zeros(size(e,1));
    elayer = zeros(Int, size(e,1));
    es = edgelengths(collect(1:size(e,1)), v, e, vertex_weight);
    elayer = vec(minimum(hcat(vlayer[e[:,1]],vlayer[e[:,2]]), 2))
    minne = 200;    maxne = 400;# 36
    
    dumparray(es, "es")
    dumparray(elayer, "elayer")
    dumparray(desiredes, "desiredes")
    dumparray(v2t, "v2t")
    dumparray(v2e, "v2e")
    # Let's get down to business
    Faile =  Int[];
    pass =1;
    while true
        availe, currvlayer = availelist(elayer, currvlayer, minne, maxne, es, nblayer, desiredes, Faile);
        print(currvlayer, ".")
        if (mod(pass,20) == 0)
            print("\n")
        end
        if (length(availe)==0)  
            break; # Done. Hallelujah!
        end 
        while true
            eL = sortedelist(elayer, es, desiredes, availe);
            Change =  false; 
            for i=1:length(eL)
                if (collapseedge!(e, es, elayer, vlayer, t, v, v2t, v2e, vertex_weight, eL[i]))
                    Change = true;  break; 
                end
                push!(Faile, eL[i]); # do I never clean up Faile?
            end
            if  (!Change)
                break;
            end
        end
        pass = pass + 1;
    end
    # # Cleanup
    print("\n")
    
    t,v,tmid =  delete_ent(t,v,tmid);
    
end

function collapseedge!(e::Array{Int, 2}, es::Vector{Float64}, elayer::Vector{Int}, vlayer::Vector{Int}, t::Array{Int, 2}, v::Array{Float64, 2}, v2t::Vector{Vector{Int}}, v2e::Vector{Vector{Int}}, vertex_weight::Vector{Float64}, dei::Int)
    result = false;
    # if the operation would result in inverted tetrahedra, cancel it
    de = e[dei, [1,2]];
    if anynegvol(t[v2t[de[2]],:], v, de[2], v[de[1],:])
        de = e[dei, [2,1]];# Try the edge the other way
        if anynegvol(t[v2t[de[2]],:], v, de[2], v[de[1],:])
            return result; # the collapse failed
        end
    end
    vi1 = de[1]; vi2 = de[2];
    mtl = unique(vcat(v2t[vi1],v2t[vi2]));
    # Modify t: switch the references to the replaced vertex vi2
    for k = 1:length(mtl)
        for i = 1:4
            if t[mtl[k], i] == vi2
                t[mtl[k], i] = vi1; break
            end 
        end
    end
    mel = unique(vcat(v2e[vi1],v2e[vi2]));
    # Modify e: switch the references to the replaced vertex vi2
    for k=1:length(mel)
        for i = 1:2
            if e[mel[k], i] == vi2
                e[mel[k], i] = vi1; break
            end 
        end
    end
    # Delete the collapsed tetrahedra from the vertex-to-tet map
    dtl = intersect(v2t[vi1],v2t[vi2]); # list of tets that connect the verts on the collapsed edge
    vl = unique(t[dtl,:]); # vertices incident on the collapsed tetrahedra
    for i=1:length(vl) # Delete the collapsed tetrahedra
        filter!(p -> !(p in dtl), v2t[vl[i]])
    end
    # Delete edges which are merged by the collapse
    del = v2e[vi2]; # vi2 is the vertex that is to be deleted
    for k =1:length(del)
        i = del[k];
        (e[i,1] == vi2) ? ov = e[i,2] : ov = e[i,1]; 
        if ov in vl
            e[i,:] = 0;# mark as deleted
            es[i] = Inf;# indicate the deleted edge
            elayer[i] = 0;# indicate the deleted edge
        end
    end
    t[dtl,:] = 0;# Mark deleted tetrahedra
    e[dei,:] = 0;# Mark deleted edge
    # Update the vertex-2-tet  map
    v2t[vi1]= setdiff(mtl, dtl);
    v2t[vi2]= [];# this vertex is gone
    # Update the vertex-2-edge  map
    v2e[vi1]= setdiff(mtl, dei);
    v2e[vi2]= [];# this vertex is gone
    #     v(vi1,:) = nv1;# new vertex location
    v[vi2,:] .= Inf;# Indicate invalid vertex
    # update edge lengths
    for k=1:length(v2e[vi1])
        i=v2e[vi1][k];
        if (e[i,1]==0)
            es[i] = Inf;# indicate the deleted edge
            elayer[i] = 0;
        else
            es[i] = elength(v[e[i,1],:], v[e[i,2],:], vertex_weight[e[i,1]], vertex_weight[e[i,2]])
            elayer[i] = minimum(vlayer[e[i,:]]);
        end
    end
    es[dei] = Inf;# indicate the deleted edge
    elayer[dei] = 0; # indicate the deleted edge
    return  true;
end

function edgelengths(ens::Vector{Int}, v::Array{Float64, 2}, e::Array{Int, 2}, vertex_weight::Vector{Float64})
    eLengths = zeros(length(ens))
    for i = 1:length(ens)
        en = ens[i]
        eLengths[i] = elength(v[e[en,1],:], v[e[en,2],:], vertex_weight[e[en,1]], vertex_weight[e[en,2]])
    end
    return eLengths
end

# Length of the edge between 2 vertices
function elength(p1::Vector{T}, p2::Vector{T}, vw1::T, vw2::T) where {T}
    p = p2 - p1;
    return max(vw1, vw2) * norm(p);
end

function sortedelist(elayer::Vector{Int}, es::Vector{Float64}, desiredes::Vector{Float64}, availe::Vector{Int})
    eList = deepcopy(availe); fill!(eList, 0)
    n=0;
    for i11=1:length(availe)
        k11=availe[i11];
        if (elayer[k11] >0)
            if  es[k11] < desiredes[elayer[k11]]
                n=n+1; eList[n] =k11;
            end
        end
    end
    return eList[1:n];
end

function availelist(elayer::Vector{Int}, currvlayer::Int, minnt::Int, maxnt::Int, es::Vector{Float64}, nblayer::Int, desiredes::Vector{Float64}, Faile::Vector{Int})
    eList= []; newcurrvlayer = 0
    for layer = currvlayer:-1:nblayer+1 # This can be changed to allow for more or less coarsening
        availe = setdiff(find(p -> layer <= p, elayer), Faile);
        eList = sortedelist(elayer, es, desiredes, availe);
        newcurrvlayer = layer;
        if (length(eList) >= minnt)
            break;
        end
    end
    availe = eList;
    availe = availe[1:min(length(availe), maxnt)];
    return availe, newcurrvlayer
end


function t4_e2(t::Array{Int, 2})
    ec = [  1  2
            2  3
            3  1
            4  1
            4  2
            4  3];
    e = vcat(t[:,ec[1,:]], t[:,ec[2,:]], t[:,ec[3,:]], t[:,ec[4,:]], t[:,ec[5,:]], t[:,ec[6,:]])
    e = sort(e, 2);
    ix = sortperm(e[:,1]);
    e = e[ix,:];
    ue = deepcopy(e)
    i = 1;
    n=1;
    while n <= size(e,1)
        c = ue[n,1];
        m = n+1;
        while m <= size(e,1)
            if (ue[m,1] != c)
                break; 
            end
            m = m+1;
        end
        us = unique(ue[n:m-1,2], 1);# the transcription below is a lot faster
        ls =length(us);
        e[i:i+ls-1,1] = c;
        e[i:i+ls-1,2] = s;
        i = i+ls;
        n = m;
    end
    e = e[1:i-1,:];
end

"""
    anynegvol(t, v, whichv, v1)

Check that the new location `v1` of the vertex `whichv` does not result 
in a negative volume  for any of the tetrahedra connected to `whichv`.  

This is a heavily used function, and hence speed and absence of
memory allocation is at a premium.
"""
function anynegvol(t::Array{Int,2}, v::Array{Float64,2}, whichv::Int, v1::Array{Float64,1})
    Volume6 = 0.0
    for iS1 = 1:size(t,1)
        i1, i2, i3, i4 = t[iS1,:]; # nodes of the tetrahedron
        if (i1 == whichv) 
            @inbounds Volume6 = let
                A1 = v[i2,1]-v1[1]; 
                A2 = v[i2,2]-v1[2]; 
                A3 = v[i2,3]-v1[3]; 
                B1 = v[i3,1]-v1[1]; 
                B2 = v[i3,2]-v1[2]; 
                B3 = v[i3,3]-v1[3]; 
                C1 = v[i4,1]-v1[1]; 
                C2 = v[i4,2]-v1[2]; 
                C3 = v[i4,3]-v1[3]; 
                ((-A3*B2+A2*B3)*C1 +  (A3*B1-A1*B3)*C2 + (-A2*B1+A1*B2)*C3)
            end
        end
        if (i2 == whichv) 
            @inbounds Volume6 = let
                A1 = v1[1]-v[i1,1]; 
                A2 = v1[2]-v[i1,2]; 
                A3 = v1[3]-v[i1,3]; 
                B1 = v[i3,1]-v[i1,1]; 
                B2 = v[i3,2]-v[i1,2]; 
                B3 = v[i3,3]-v[i1,3]; 
                C1 = v[i4,1]-v[i1,1]; 
                C2 = v[i4,2]-v[i1,2]; 
                C3 = v[i4,3]-v[i1,3]; 
                ((-A3*B2+A2*B3)*C1 +  (A3*B1-A1*B3)*C2 + (-A2*B1+A1*B2)*C3)
            end
        end
        if (i3 == whichv) 
            @inbounds Volume6 = let
                A1 = v[i2,1]-v[i1,1]; 
                A2 = v[i2,2]-v[i1,2]; 
                A3 = v[i2,3]-v[i1,3]; 
                B1 = v1[1]-v[i1,1]; 
                B2 = v1[2]-v[i1,2]; 
                B3 = v1[3]-v[i1,3]; 
                C1 = v[i4,1]-v[i1,1]; 
                C2 = v[i4,2]-v[i1,2]; 
                C3 = v[i4,3]-v[i1,3]; 
                ((-A3*B2+A2*B3)*C1 +  (A3*B1-A1*B3)*C2 + (-A2*B1+A1*B2)*C3)
            end
        end
        if (i4 == whichv) 
            @inbounds Volume6 = let
                A1 = v[i2,1]-v[i1,1]; 
                A2 = v[i2,2]-v[i1,2]; 
                A3 = v[i2,3]-v[i1,3]; 
                B1 = v[i3,1]-v[i1,1]; 
                B2 = v[i3,2]-v[i1,2]; 
                B3 = v[i3,3]-v[i1,3]; 
                C1 = v1[1]-v[i1,1]; 
                C2 = v1[2]-v[i1,2]; 
                C3 = v1[3]-v[i1,3]; 
                ((-A3*B2+A2*B3)*C1 +  (A3*B1-A1*B3)*C2 + (-A2*B1+A1*B2)*C3)
            end
        end
        if (Volume6 < 0)
             return true;
        end
    end
    return false;
end

function delete_ent(t::Array{Int,2}, v::Array{Float64,2}, tmid::Array{Int,1})
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
            if (!isempty(find(p -> p == 0, nn[t[i,:]])))
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