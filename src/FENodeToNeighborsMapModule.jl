"""
    FENodeToNeighborsMapModule  

Module to construct a map from finite element nodes to nodes connected to them.
"""
module FENodeToNeighborsMapModule

__precompile__(true)

using ..FENodeToFEMapModule: FENodeToFEMap

function __collect_unique_node_neighbors(ellist, conn, npe)
    totn = length(ellist) * npe
    nodes = fill(zero(eltype(conn[1])), totn)
    p = 1
    @inbounds for i in ellist
        for k in conn[i]
            nodes[p] = k 
            p += 1
        end
    end
    sort!(nodes)
    unique!(nodes)
    return nodes
end

function _unique_nodes(n2e, conn)
    npe = length(conn[1])
    empt = eltype(n2e.map[1])[]
    unique_nodes = fill(empt, length(n2e.map))
    Base.Threads.@threads for i in 1:length(n2e.map) # run this in PARALLEL
        unique_nodes[i] = __collect_unique_node_neighbors(n2e.map[i], conn, npe)
    end
    return unique_nodes
end

"""
    FENodeToNeighborsMap

Map from finite element nodes to the nodes that are connected to them by finite
elements.

!!! note

    Self references are included (each node is connected to itself).
"""
struct FENodeToNeighborsMap{IT}
    # Map as a vector of vectors.
    map::Vector{Vector{IT}}
end

"""
    FENodeToNeighborsMap(conn::Vector{NTuple{N, IT}}, nmax::FInt) where {N, IT<:Integer}

Map from finite element nodes to the finite elements connecting them.

- `conns` = connectivities as a vector of tuples
- `nmax` = largest possible node number

Example:
```
m = FENodeToNeighborsMap(fes.conn, count(fens))
```
"""
function FENodeToNeighborsMap(n2e::N2EMAP, conn::Vector{NTuple{N, IT}}) where {N2EMAP<:FENodeToFEMap, N, IT<:Integer}
    return FENodeToNeighborsMap(_unique_nodes(n2e, conn))
end

end
