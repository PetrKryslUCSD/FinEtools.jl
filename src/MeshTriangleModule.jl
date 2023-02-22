"""
    MeshTriangleModule

Module  for generation of meshes composed of triangles.
"""
module MeshTriangleModule

__precompile__(true)

import LinearAlgebra: norm
using ..FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import ..FESetModule: AbstractFESet, FESetT3, FESetT6, FESetQ4, connasarray, subset
import ..FENodeSetModule: FENodeSet
import ..MeshUtilModule: makecontainer, addhyperface!, findhyperface!, linearspace
import ..MeshModificationModule:
    meshboundary, connectednodes, fusenodes, updateconn!, compactnodes, renumberconn!
import ..MeshSelectionModule: selectelem, findunconnnodes
import ..MeshQuadrilateralModule: Q4block
import Statistics: mean

"""
    T3blockx(xs::FFltVec, ys::FFltVec, orientation::Symbol=:a)

T3 Mesh of a rectangle.
"""
function T3blockx(xs::FFltVec, ys::FFltVec, orientation::Symbol = :a)
    if (orientation != :a) && (orientation != :b)
        error("Cannot handle orientation : $orientation")
    end
    nL = length(xs) - 1
    nW = length(ys) - 1
    nnodes = (nL + 1) * (nW + 1)
    ncells = 2 * (nL) * (nW)
    xys = zeros(FFlt, nnodes, 2)
    conns = zeros(FInt, ncells, 3)
    f = 1
    for j in 1:(nW+1)
        for i in 1:(nL+1)
            xys[f, 1] = xs[i]
            xys[f, 2] = ys[j]
            f = f + 1
        end
    end
    fens = FENodeSet(xys)

    gc = 1
    for i in 1:nL
        for j in 1:nW
            f = (j - 1) * (nL + 1) + i
            if (orientation == :a)
                nn = [f (f + 1) f + (nL + 1)]
            elseif (orientation == :b)
                nn = [f (f + 1) f + (nL + 1) + 1]
            end
            conns[gc, :] = nn
            gc = gc + 1
            if (orientation == :a)
                nn = [(f + 1) f + (nL + 1) + 1 f + (nL + 1)]
            elseif (orientation == :b)
                nn = [f f + (nL + 1) + 1 f + (nL + 1)]
            end
            conns[gc, :] = nn
            gc = gc + 1
        end
    end
    fes = FESetT3(conns)
    return fens, fes
end

function T3blockx(xs::FFltMat, ys::FFltMat, orientation::Symbol = :a)
    return T3blockx(vec(xs), vec(ys), orientation)
end

"""
    T3block(Length::FFlt, Width::FFlt, nL::FInt, nW::FInt, orientation::Symbol=:a)

T3 Mesh of a rectangle.
"""
function T3block(Length::FFlt, Width::FFlt, nL::FInt, nW::FInt, orientation::Symbol = :a)
    return T3blockx(
        FFltVec(linearspace(0.0, Length, nL + 1)),
        FFltVec(linearspace(0.0, Width, nW + 1)),
        orientation,
    )
end

"""
    T3toT6(fens::FENodeSet, fes::FESetT3)

Convert a mesh of triangle T3 (three-node) to triangle T6.
"""
function T3toT6(fens::FENodeSet, fes::FESetT3)
    nedges = 3
    ec = [1 2; 2 3; 3 1]
    conns = connasarray(fes)
    # Additional node numbers are numbered from here
    newn = count(fens) + 1
    # make a search structure for edges
    edges = makecontainer()
    for i in axes(conns, 1)
        conn = conns[i, :]
        for J in 1:nedges
            ev = conn[ec[J, :]]
            newn = addhyperface!(edges, ev, newn)
        end
    end
    xyz1 = fens.xyz             # Pre-existing nodes
    # Allocate for vertex nodes plus edge nodes plus face nodes
    xyz = zeros(FFlt, newn - 1, size(fens.xyz, 2))
    xyz[1:size(xyz1, 1), :] = xyz1 # existing nodes are copied over
    # calculate the locations of the new nodes
    # and construct the new nodes
    for i in keys(edges)
        C = edges[i]
        for J in eachindex(C)
            ix = vec([item for item in C[J].o])
            push!(ix, i)
            xyz[C[J].n, :] = mean(xyz[ix, :], dims = 1)
        end
    end
    # construct new geometry cells
    nconns = zeros(FInt, size(conns, 1), 6)
    nc = 1
    for i in axes(conns, 1)
        conn = conns[i, :]
        econn = zeros(FInt, 1, nedges)
        for J in 1:nedges
            ev = conn[ec[J, :]]
            h, n = findhyperface!(edges, ev)
            econn[J] = n
        end
        nconns[nc, :] = vcat(vec(conn), vec(econn))
        nc = nc + 1
    end
    fens = FENodeSet(xyz)
    fes = FESetT6(nconns)
    return fens, fes
end

"""
    T6block(Length::FFlt, Width::FFlt, nL::FInt, nW::FInt, orientation::Symbol=:a)

Mesh of a rectangle of T6 elements.
"""
function T6block(Length::FFlt, Width::FFlt, nL::FInt, nW::FInt, orientation::Symbol = :a)
    fens, fes = T3block(Length, Width, nL, nW, orientation)
    fens, fes = T3toT6(fens, fes)
end

"""
    T6blockx(xs::FFltVec, ys::FFltVec, orientation::Symbol=:a)

Graded mesh of a 2-D block of T6 finite elements.
"""
function T6blockx(xs::FFltVec, ys::FFltVec, orientation::Symbol = :a)
    fens, fes = T3blockx(xs, ys, orientation)
    fens, fes = T3toT6(fens, fes)
end

"""
    Q4toT3(fens::FENodeSet, fes::FESetQ4, orientation::Symbol=:default)

Convert a mesh of quadrilateral Q4's to two T3 triangles  each.
"""
function Q4toT3(fens::FENodeSet, fes::FESetQ4, orientation::Symbol = :default)
    numbering1 = ([1, 2, 3], [1, 3, 4])
    numbering2 = ([1, 2, 4], [3, 4, 2])
    numbering = numbering1
    if orientation == :alternate
        numbering = numbering2
    end
    nedges = 4
    nconns = zeros(FInt, 2 * count(fes), 3)
    conns = connasarray(fes)
    if orientation == :random
        draw = rand(Bool, size(conns, 1))
        nc = 1
        for i in axes(conns, 1)
            conn = conns[i, :]
            if draw[i]
                nconns[nc, :] = conn[numbering1[1]]
                nc = nc + 1
                nconns[nc, :] = conn[numbering1[2]]
                nc = nc + 1
            else
                nconns[nc, :] = conn[numbering2[1]]
                nc = nc + 1
                nconns[nc, :] = conn[numbering2[2]]
                nc = nc + 1
            end
        end
    else
        nc = 1
        for i in axes(conns, 1)
            conn = conns[i, :]
            nconns[nc, :] = conn[numbering[1]]
            nc = nc + 1
            nconns[nc, :] = conn[numbering[2]]
            nc = nc + 1
        end
    end
    nfes = FESetT3(nconns)
    return fens, nfes            # I think I should not be overwriting the input!
end

"""
    T3refine(fens::FENodeSet,fes::FESetT3)

Refine a mesh of 3-node triangles by quadrisection.
"""
function T3refine(fens::FENodeSet, fes::FESetT3)
    fens, fes = T3toT6(fens, fes)
    nconn = zeros(FInt, 4 * size(fes.conn, 1), 3)
    conns = connasarray(fes)
    nc = 1
    for i in axes(conns, 1)
        c = conns[i, :]
        nconn[nc, :] = c[[1, 4, 6]]
        nc = nc + 1
        nconn[nc, :] = c[[2, 5, 4]]
        nc = nc + 1
        nconn[nc, :] = c[[3, 6, 5]]
        nc = nc + 1
        nconn[nc, :] = c[[4, 5, 6]]
        nc = nc + 1
    end
    nfes = FESetT3(nconn)
    return fens, nfes            # I think I should not be overwriting the input!
end

"""
    T3annulus(rin::FFlt, rex::FFlt, nr::FInt, nc::FInt, Angl::FFlt, orientation::Symbol=:a)

Mesh of an annulus segment.

Mesh of an annulus segment, centered at the origin, with internal radius `rin`,
and  external radius `rex`, and  development angle `Angl` (in radians). Divided
into elements: nr, nc in the radial and circumferential direction respectively.
"""
function T3annulus(
    rin::FFlt,
    rex::FFlt,
    nr::FInt,
    nc::FInt,
    Angl::FFlt,
    orientation::Symbol = :a,
)
    trin = min(rin, rex)
    trex = max(rin, rex)
    fens, fes = T3block(trex - trin, Angl, nr, nc, orientation)
    xy = fens.xyz
    for i in 1:count(fens)
        r = trin + xy[i, 1]
        a = xy[i, 2]
        xy[i, :] = [r * cos(a) r * sin(a)]
    end
    fens.xyz = xy
    return fens, fes
end

"""
    T6annulus(rin::FFlt, rex::FFlt, nr::FInt, nc::FInt, Angl::FFlt, orientation::Symbol=:a)

Mesh of an annulus segment.

Mesh of an annulus segment, centered at the origin, with internal radius `rin`,
and  external radius `rex`, and  development angle `Angl` (in radians). Divided
into elements: nr, nc in the radial and circumferential direction respectively.
"""
function T6annulus(
    rin::FFlt,
    rex::FFlt,
    nr::FInt,
    nc::FInt,
    Angl::FFlt,
    orientation::Symbol = :a,
)
    trin = min(rin, rex)
    trex = max(rin, rex)
    fens, fes = T3block(trex - trin, Angl, nr, nc, orientation)
    fens, fes = T3toT6(fens, fes)
    xy = fens.xyz
    for i in 1:count(fens)
        r = trin + xy[i, 1]
        a = xy[i, 2]
        xy[i, :] = [r * cos(a) r * sin(a)]
    end
    fens.xyz = xy
    return fens, fes
end


"""
    T3circlen(radius::FFlt, nperradius)

Mesh of a quarter circle with a given number of elements per radius.

The parameter `nperradius` should be an even 
number; if that isn't so is adjusted to by adding one. 
"""
function T3circlen(radius::FFlt, nperradius)
    fens, fes = T3block(1.0, 1.0, 1, 1)
    fens.xyz[1, 1] = 1.0
    fens.xyz[1, 2] = 0.0
    fens.xyz[2, 1] = sqrt(2.0) / 2
    fens.xyz[2, 2] = sqrt(2.0) / 2
    fens.xyz[3, 1] = 0.0
    fens.xyz[3, 2] = 0.0
    fens.xyz[4, 1] = 0.0
    fens.xyz[4, 2] = 1.0
    tolerance = 1.0 / nperradius / 10
    n = 1
    while n < nperradius
        fens, fes = T3refine(fens, fes)
        bfes = meshboundary(fes)
        lx = selectelem(fens, bfes; box = Float64[0 0 -Inf Inf], inflate = tolerance)
        ly = selectelem(fens, bfes; box = Float64[-Inf Inf 0 0], inflate = tolerance)
        lc = setdiff(1:count(bfes), vcat(lx, ly))
        l1 = connectednodes(subset(bfes, lc))
        for j in l1
            d = norm(fens.xyz[j, :])
            fens.xyz[j, :] .*= 1.0 / d
        end
        n = n * 2
    end
    fens.xyz .*= radius
    return fens, fes
end

"""
    T3circleseg(angle::FFlt, radius::FFlt, nperradius, ncircumferentially, orientation::Symbol=:a)

Mesh of a segment of a circle.

The subtended angle is `angle` in radians. The orientation: refer to `T3block`.
"""
function T3circleseg(
    angle::FFlt,
    radius::FFlt,
    ncircumferentially,
    nperradius,
    orientation::Symbol = :a,
)
    fens, fes = T3block(angle, radius, ncircumferentially, nperradius, orientation)
    for i in 1:count(fens)
        a = angle - fens.xyz[i, 1]
        r = fens.xyz[i, 2]
        fens.xyz[i, :] .= (r * cos(a), r * sin(a))
    end
    tolerance = radius / max(nperradius, ncircumferentially) / 100
    fens, newn = fusenodes(fens, fens, tolerance)
    updateconn!(fes, newn)
    connected = findunconnnodes(fens, fes)
    fens, newn = compactnodes(fens, connected)
    fes = renumberconn!(fes, newn)
    keep = fill(true, count(fes))
    ca = connasarray(fes)
    for j in 1:count(fes)
        keep[j] = (length(unique(ca[j, :])) == 3)
    end
    l = findall(x -> x == true, keep)
    fes = subset(fes, l)
    return fens, fes
end


end
