module MortarModule

__precompile__(true)

using LinearAlgebra
using SparseArrays

using ..FENodeSetModule: FENodeSet
using ..FESetModule: count
using ..IntegRuleModule: TriRule, GaussRule
using ..IntegDomainModule: IntegDomain
using ..NodalFieldModule: NodalField
using ..ElementalFieldModule: ElementalField
using ..FEMMBaseModule: FEMMBase, bilform_masslike, bilform_dot
using ..FESetModule: FESetT3, FESetT6
using ..DataCacheModule: DataCache
using ..FieldModule: numberdofs!

export common_refinement

function clip_polygon(subject::Vector{Vector{Float64}},
                      clip::Vector{Vector{Float64}}; eps=1e-12)

    isempty(subject) && return Vector{Vector{Float64}}[]
    isempty(clip)    && return Vector{Vector{Float64}}[]

    same(p, q) = norm(p - q) <= eps

    function poly_normal(poly)
        n = zeros(3)
        m = length(poly)
        for i in 1:m
            p = poly[i]
            q = poly[mod1(i + 1, m)]
            n += cross(p, q)
        end
        nn = norm(n)
        nn <= eps && error("Degenerate clip polygon")
        return n / nn
    end

    # Orthogonal projection of point P onto the plane of clip polygon.
    function project_to_clip_plane(P, P0, n)
        return P - dot(P - P0, n) * n
    end

    function side(P, A, B, n)
        return dot(cross(B - A, P - A), n)
    end

    function intersect_seg_clipedge(S, E, A, B, n)
        sS = side(S, A, B, n)
        sE = side(E, A, B, n)
        den = sS - sE

        if abs(den) <= eps
            return nothing
        end

        t = sS / den
        return S + t * (E - S)
    end

    function push_unique!(out, p)
        for q in out
            same(p, q) && return
        end
        push!(out, p)
    end

    function clean(poly)
        out = Vector{Vector{Float64}}()

        for p in poly
            if isempty(out) || !same(out[end], p)
                push!(out, p)
            end
        end

        if length(out) > 1 && same(out[1], out[end])
            pop!(out)
        end

        out2 = Vector{Vector{Float64}}()
        for p in out
            push_unique!(out2, p)
        end

        return out2
    end

    n = poly_normal(clip)
    P0 = clip[1]

    # This is the key change:
    # project the subject polygon onto the clip polygon plane.
    output = [
        project_to_clip_plane(p, P0, n)
        for p in subject
    ]

    output = clean(output)
    length(output) < 3 && return Vector{Vector{Float64}}[]

    for i in 1:length(clip)
        A = clip[i]
        B = clip[mod1(i + 1, length(clip))]

        input = output
        output = Vector{Vector{Float64}}()
        isempty(input) && break

        S = input[end]
        sS = side(S, A, B, n)
        Sin = sS >= -eps

        for E in input
            sE = side(E, A, B, n)
            Ein = sE >= -eps

            if Ein
                if !Sin
                    X = intersect_seg_clipedge(S, E, A, B, n)
                    X !== nothing && push!(output, X)
                end
                push!(output, E)
            elseif Sin
                X = intersect_seg_clipedge(S, E, A, B, n)
                X !== nothing && push!(output, X)
            end

            S = E
            sS = sE
            Sin = Ein
        end

        output = clean(output)
    end

    output = clean(output)

    # Numerical safety: force final vertices exactly back onto B/clip plane.
    output = [
        project_to_clip_plane(p, P0, n)
        for p in output
    ]

    output = clean(output)

    if length(output) < 3
        return Vector{Vector{Float64}}[]
    end

    return output
end

# ##############################################################################
function get_node_id(x::Vector{Float64}, node_map, XU,
                     ai, ax, IA, JA, VA, 
                     bi, bx, IB, JB, VB; order=1, dim_u=1)
    # TODO: cannot use anything in the negative

    return get!(node_map, abs.(round.(x; digits=5))) do
        # create new node
        push!(XU, x)
        
        baryA = barycentre(x, ax)
        # 
        nids_a = vcat([[dim_u*(size(XU,1)-1) + j for j in 1:dim_u] for k in 1:length(ai)]...)
        nids_b = vcat([[dim_u*(size(XU,1)-1) + j for j in 1:dim_u] for k in 1:length(bi)]...)

        dofs_a = vcat([[dim_u*(ai[k]-1) + j for j in 1:dim_u] for k in 1:length(ai)]...)
        baryAs = [baryA[k] for k in 1:length(ai) for j in 1:dim_u]

        push!(IA, nids_a...)
        push!(JA, dofs_a...)
        push!(VA, baryAs...)
        if order == 1
            baryB = barycentre(x, bx)
            baryBs = [baryB[k] for k in 1:length(bi) for j in 1:dim_u]
            dofs_b = vcat([[dim_u*(bi[k]-1) + j for j in 1:dim_u] for k in 1:length(bi)]...)
            push!(IB, nids_b...)
            push!(JB, dofs_b...)
            push!(VB, baryBs...)
        end
        size(XU,1)

    end
end
# ###############################################################################
struct Grid
    origin::Vector{Float64}
    h::Float64
    cells::Dict{NTuple{3,Int}, Vector{Int}}
end

cell_index(g::Grid, x::Vector{Float64}) = (
    floor(Int, (x[1]-g.origin[1])/g.h),
    floor(Int, (x[2]-g.origin[2])/g.h),
    floor(Int, (x[3]-g.origin[3])/g.h)
)

function aabb(X::Matrix{Float64}, poly_conn::Vector{Int})
    pts = X[poly_conn, :]
    mn = vec(minimum(pts, dims=1))
    mx = vec(maximum(pts, dims=1))
    return mn, mx
end

function build_grid(X::Matrix{Float64}, Conn::Matrix{Int}; h::Float64, pad::Float64=0.0)
    mn = vec(minimum(X, dims=1)) .- pad
    g = Grid(mn, h, Dict{NTuple{3,Int}, Vector{Int}}())
    for k in 1:size(Conn,1)
        tri = Conn[k, :]
        aabb_mn, aabb_mx = aabb(X, tri)
        aabb_mn .-= pad
        aabb_mx .+= pad
        i0 = cell_index(g, aabb_mn)
        i1 = cell_index(g, aabb_mx)
        for i in i0[1]:i1[1], j in i0[2]:i1[2], l in i0[3]:i1[3]
            key = (i,j,l)
            push!(get!(g.cells, key, Int[]), k)
        end
    end
    return g
end

function query_grid(g::Grid, aabb_mn::Vector{Float64}, aabb_mx::Vector{Float64})
    i0 = cell_index(g, aabb_mn)
    i1 = cell_index(g, aabb_mx)

    cand = Int[]
    for i in i0[1]:i1[1], j in i0[2]:i1[2], l in i0[3]:i1[3]
        key = (i,j,l)
        if haskey(g.cells, key)
            append!(cand, g.cells[key])
        end
    end
    sort!(cand)
    unique!(cand)
    return cand
end
# ###############################################################################

function barycentre(x, xas::Matrix)
    if size(xas, 1) == 3
        a = xas[1, :]
        b = xas[2, :]
        c = xas[3, :]

        v0 = b - a
        v1 = c - a

        d00 = dot(v0, v0)
        d01 = dot(v0, v1)
        d11 = dot(v1, v1)

        denom = d00 * d11 - d01 * d01
        @assert denom != 0 "Triangle is degenerate"
        v2 = x - a

        d20 = dot(v2, v0)
        d21 = dot(v2, v1)

        lam2 = (d11 * d20 - d01 * d21) / denom
        lam3 = (d00 * d21 - d01 * d20) / denom
        lam1 = 1.0 - lam2 - lam3

        bary = [lam1, lam2, lam3]

        return bary
    elseif size(xas, 1) == 4
        a = xas[1, :]
        b = xas[2, :]
        c = xas[3, :]
        d = xas[4, :]

        # Bilinear map:
        # X(s,t) = (1-s)(1-t)a + s(1-t)b + st c + (1-s)t d
        #
        # Equivalently:
        # X(s,t) = a + s(b-a) + t(d-a) + s*t(c-b-d+a)

        e1 = b - a
        e2 = d - a
        e3 = c - b - d + a

        # good default for clipped points inside the element
        s = 0.5
        t = 0.5
        for iter in 1:20
            r = a + s*e1 + t*e2 + (s*t)*e3 - x

            Js = e1 + t*e3
            Jt = e2 + s*e3

            J = hcat(Js, Jt)   # 3 x 2
            delta = J \ r          # least-squares, no singular planar 3x3 solve

            s -= delta[1]
            t -= delta[2]

            if norm(delta) < 1e-12 || norm(r) < 1e-12
                break
            end
        end

        lam1 = (1 - s) * (1 - t)
        lam2 = s * (1 - t)
        lam3 = s * t
        lam4 = (1 - s) * t

        return [lam1, lam2, lam3, lam4]
    else
        error("Unsupported element with $(size(xas,1)) nodes")
    end
end 
# ###############################################################################
function common_refinement(fensA, fesA,
                            fensB, fesB; 
                            h = 0.1, lam_order = 1, tri_order = 1,
                             triangulation_type = "naive" , dim_u=1)

    XA = fensA.xyz
    connA = stack(fesA.conn, dims=1)
    XB = fensB.xyz
    connB = stack(fesB.conn, dims=1)
    if size(fensA.xyz, 2) !=3
        XA = hcat(XA, zeros(size(XA, 1)))
    end
    if size(fensB.xyz, 2) !=3
        XB = hcat(XB, zeros(size(XB, 1)))
    end

    nA = size(connA,1)
    nB = size(connB,1)
    pad = 1e-2
    # 
    gridB = build_grid(XB, connB; h=h, pad=pad)

    parentA = Int[]
    parentB = Int[]

    # containers for barycentric points
    IA = Int[]; JA = Int[]; VA = Float64[]
    IB = Int[]; JB = Int[]; VB = Float64[]

    SIA = Int[]; SJA = Int[]; SVA = Float64[]
    SIB = Int[]; SJB = Int[]; SVB = Float64[]


    XU = Vector{Vector{Float64}}()
    connU = Array{Int}[]
    node_map = Dict{Vector{Float64},Int}()
    for i in 1:nA
        count = 0
        ai = connA[i,:]
        aabb_mn, aabb_mx = aabb(XA, ai)
        # pad = 1e-6
        aabb_mn .-= pad
        aabb_mx .+= pad
        cands = query_grid(gridB, aabb_mn, aabb_mx)
        # println("aabb_min: $aabb_mn, aabb_mx: $aabb_mx")
        # @info "Cands for element $i in A: $(length(cands))"
        for j in cands
            count += 1
           
            
            bi = connB[j,:]
            ai = connA[i,:]
            ax = XA[connA[i,:], :]
            bx = XB[connB[j,:], :]



            # Clipping: change to custom implementation is needed
            # clipped = clip(PolyArea([(XA[connA[i, k], 1], XA[connA[i, k], 2]) for k in 1:size(connA,2)]), 
            #     PolyArea([(XB[connB[j, k], 1], XB[connB[j, k], 2]) for k in 1:size(connB,2)]), 
            #     SutherlandHodgmanClipping())
            
            clipped = clip_polygon([vec(XA[connA[i, k], :]) for k in 1:size(connA,2)],
                                  [vec(XB[connB[j, k], :]) for k in 1:size(connB,2)])
                             
            if isnothing(clipped)
                continue
            end
            

            # nv = length(clipped.rings[1].vertices)
            nv = length(clipped)
            conn = Vector{Int}()
            for k in 1:nv
                # vs = [clipped.rings[1].vertices[k].coords.x.val, clipped.rings[1].vertices[k].coords.y.val, 0.0]
                
                vs = [clipped[k][1], clipped[k][2], clipped[k][3]]
                push!(conn, get_node_id(vs, node_map, XU,
                                        ai, ax, IA, JA, VA, 
                                        bi, bx, IB, JB, VB; order=lam_order, dim_u=dim_u))
            end
            conn = unique(conn)
            nv = length(conn)
            # # Triangulation
            if triangulation_type=="naive"
                if length(conn)>=3
                    for k in 3:nv
                        conn_cuurent = [conn[1], conn[k-1], conn[k]]

                        if tri_order == 2
                            # add midpoints of edges
                            mid12 = (XU[conn[1]] + XU[conn[k-1]]) / 2
                            mid23 = (XU[conn[k-1]] + XU[conn[k]]) / 2
                            mid31 = (XU[conn[k]] + XU[conn[1]]) / 2

                            mid12_id = get_node_id(mid12, node_map, XU,
                                                    ai, ax, IA, JA, VA, 
                                                    bi, bx, IB, JB, VB; order=lam_order, dim_u=dim_u)
                            mid23_id = get_node_id(mid23, node_map, XU,
                                                    ai, ax, IA, JA, VA, 
                                                    bi, bx, IB, JB, VB; order=lam_order, dim_u=dim_u)
                            mid31_id = get_node_id(mid31, node_map, XU,
                                                    ai, ax, IA, JA, VA, 
                                                    bi, bx, IB, JB, VB; order=lam_order, dim_u=dim_u)

                            conn_cuurent = [conn[1], conn[k-1], conn[k], mid12_id, mid23_id, mid31_id]
                        end
                        push!(connU, conn_cuurent )
                        push!(parentA, i)
                        push!(parentB, j)
                        if lam_order==0
                            push!(IB, [dim_u*(size(connU,1)-1) + ind for ind in 1:dim_u]...)
                            push!(JB, [dim_u*(j-1) + ind for ind in 1:dim_u]...)
                            push!(VB, [ 1.0 for j in 1:dim_u]...)
                        end
                    end
                end
            elseif triangulation_type=="cp"
                # get midpoint of clipped polygon for use as Steiner point in triangulation
                centroid = zeros(3)
                for k in 1:nv
                    centroid .+= XU[conn[k]]
                end
                centroid ./= nv
                cnid = get_node_id(centroid, node_map, XU,
                            ai, ax, IA, JA, VA, 
                            bi, bx, IB, JB, VB; order=lam_order, dim_u=dim_u)
                # triangulate using centroid as Steiner 
                for k in 1:nv
                    conn_cuurent = [conn[k], conn[mod1(k+1, nv)], cnid]

                    if tri_order == 2
                        # add midpoints of edges
                        mid12 = (XU[conn[k]] + XU[conn[mod1(k+1, nv)]]) / 2
                        mid23 = (XU[conn[mod1(k+1, nv)]] + XU[cnid]) / 2
                        mid31 = (XU[cnid] + XU[conn[k]]) / 2

                        mid12_id = get_node_id(mid12, node_map, XU,
                                                ai, ax, IA, JA, VA, 
                                                bi, bx, IB, JB, VB; order=lam_order, dim_u=dim_u)
                        mid23_id = get_node_id(mid23, node_map, XU,
                                                ai, ax, IA, JA, VA, 
                                                bi, bx, IB, JB, VB; order=lam_order, dim_u=dim_u)
                        mid31_id = get_node_id(mid31, node_map, XU,
                                                ai, ax, IA, JA, VA, 
                                                bi, bx, IB, JB, VB; order=lam_order, dim_u=dim_u)

                        conn_cuurent = [conn[k], conn[mod1(k+1,nv)], cnid , mid12_id, mid23_id, mid31_id]
                    end
                    push!(connU, conn_cuurent )
                    push!(parentA, i)
                    push!(parentB, j)
                    if lam_order==0
                        push!(IB, [dim_u*(size(connU,1)-1) + ind for ind in 1:dim_u]...)
                        push!(JB, [dim_u*(j-1) + ind for ind in 1:dim_u]...)
                        push!(VB, [ 1.0 for j in 1:dim_u]...)
                    end
                end
            end
        end
    end
    # making dimensions consistent for C/D
    

    nXu = length(XU)
    nUe = length(connU)

    nA_dofs = dim_u * size(XA, 1)
    nB_dofs = lam_order == 0 ? dim_u * size(connB, 1) : dim_u * size(XB, 1)

    PiA = sparse(IA, JA, VA, dim_u * nXu, nA_dofs)

    if lam_order == 0
        PiB = sparse(IB, JB, VB, dim_u * nUe, dim_u * size(connB, 1))
    else
        PiB = sparse(IB, JB, VB, dim_u * nXu, dim_u * size(XB, 1))
    end

    # PiA = 0
    # PiB = 0


    # 
    if tri_order ==1
        fesu = FESetT3(stack(connU, dims=1))
    elseif tri_order == 2
        fesu = FESetT6(stack(connU, dims=1))
    end
    fensu = FENodeSet(stack(XU, dims=1))
    # if tri_order==2
    #     fensu, fesu = T3toT6(fensu, fesu)
    # end
    geomu = NodalField(fensu.xyz)
    uu = NodalField(zeros(size(fensu.xyz, 1), dim_u))
    numberdofs!(uu)
    femmu = FEMMBase(IntegDomain(fesu, TriRule(3)))

    if lam_order ==1
        # @infiltrate
        M = bilform_dot(femmu, geomu, uu, DataCache(LinearAlgebra.I(dim_u)))

        # onevec = ones(size(M, 1))
    elseif lam_order == 0
        M = bilform_masslike(femmu, geomu, uu, DataCache(LinearAlgebra.I(dim_u)))
    end
    C = PiB' * M * PiA

    

    meta = Dict(
        "XA" => XA,
        "connA" => connA,
        "XB" => XB,
        "connB" => connB,
        "fes_u" => fesu,
        "fens_u" => fensu,
        "parentA" => parentA,
        "parentB" => parentB,
        "PiA" => PiA,
        "PiB" => PiB,
        "M" => M,
    )

    return C, meta
end


end