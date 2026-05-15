"""
    MeshLineModule

Module  for generation of meshes composed of line (curve) elements.
"""
module MeshLineModule

__precompile__(true)

import ..FENodeSetModule: FENodeSet
import ..FESetModule: FESetL2, FESetL3
import ..MeshUtilModule: linearspace
import Statistics: mean

"""
    L2block(Length::T, nL::IT) where {T<:Number, IT<:Integer}

Mesh of a 1-D block of L2 finite elements.
"""
function L2block(Length::T, nL::IT) where {T<:Number,IT<:Integer}
    fens, fes = L2blockx(collect(dropdims(linearspace(0.0, Length, nL + 1)', dims = 1)))
end

"""
    L2blockx(xs::Vector{T}) where {T<:Number}

Graded mesh of a 1-D block, L2 finite elements.
"""
function L2blockx(xs::Vector{T}) where {T<:Number}
    xyz = reshape(xs, length(xs), 1)
    ncells = length(xs) - 1

    # create the nodes
    fens = FENodeSet(xyz)
    # Create the finite elements
    fes = FESetL2([(1:ncells) (2:(ncells+1))])

    return fens, fes
end

"""
    L2blockx(xs::Vector{T}, ys::Vector{T}) where {T<:Number}

Graded mesh of a 1-D block, L2 finite elements.

Embedded in two dimensions.
"""
function L2blockx(xs::Vector{T}, ys::Vector{T}) where {T<:Number}
    @assert length(xs) == length(ys) "The vectors of coordinates must have the same length."
    xyz = reshape([xs; ys], length(xs), 2)
    ncells = length(xs) - 1

    # create the nodes
    fens = FENodeSet(xyz)
    # Create the finite elements
    fes = FESetL2([(1:ncells) (2:(ncells+1))])

    return fens, fes
end

"""
    L2blockx(xs::Vector{T}, ys::Vector{T}, zs::Vector{T}) where {T<:Number}

Graded mesh of a 1-D block, L2 finite elements.

Embedded in three dimensions.
"""
function L2blockx(xs::Vector{T}, ys::Vector{T}, zs::Vector{T}) where {T<:Number}
    @assert length(xs) == length(ys) "The vectors of coordinates must have the same length."
    @assert length(xs) == length(zs) "The vectors of coordinates must have the same length."
    xyz = reshape([xs; ys; zs], length(xs), 3)
    ncells = length(xs) - 1

    # create the nodes
    fens = FENodeSet(xyz)
    # Create the finite elements
    fes = FESetL2([(1:ncells) (2:(ncells+1))])

    return fens, fes
end

"""
    L3blockx(xs::Vector{T}) where {T<:Number}

Graded mesh of a 1-D block, L2 finite elements.
"""
function L3blockx(xs::Vector{T}) where {T<:Number}
    fens, fes = L2blockx(xs)
    nxyz = zeros(count(fes), size(fens.xyz, 2))
    nconn = zeros(Int, count(fes), 3)
    N = count(fens)
    for i in eachindex(fes)
        N = N + 1
        nxyz[i, :] = mean(fens.xyz[[k for k in fes.conn[i]], :], dims = 1)
        nconn[i, :] = vcat([k for k in fes.conn[i]], [N])
    end
    fens = FENodeSet([fens.xyz; nxyz])
    # Create the finite elements
    fes = FESetL3(nconn)

    return fens, fes
end

end
