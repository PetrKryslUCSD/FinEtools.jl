"""
    FENodeToFEMapModule

Module to construct a map from finite element nodes to the finite elements.
"""
module FENodeToFEMapModule

using ..FESetModule: AbstractFESet

__precompile__(true)

"""
    FENodeToFEMap

Map from finite element nodes to the finite elements connecting them.

For each  node referenced in the connectivity of the finite element set on
input, the numbers of the individual finite elements that reference that node is
stored in an array in the array `map`.

        Example: 
```
fes.conn= [7,6,5;
            4,1,3;
            3,7,5];
The map reads
    map[1] = [2];
    map[2] = [];#  note that node number 2 is not referenced by the connectivity
    map[3] = [2,3];
    map[4] = [2];
    map[5] = [1,3];
    map[6] = [1];
    map[7] = [1,3];
```
The individual elements from the connectivity that reference node number 5 are 1
and 3, so that `fes.conn(map[5],:) `includes all the nodes that are connected to
node 5 (including node 5 itself).
"""
struct FENodeToFEMap{IT}
    # Map as a vector of vectors.
    map::Vector{Vector{IT}}
end

function _makemap(conn, fesindexes, nmax::IT) where {IT<:Integer}
    map = Array{Vector{IT}}(undef, nmax)
    @inbounds for i in eachindex(map)
        map[i] = IT[]  # initially empty arrays
    end
    @inbounds for j in fesindexes
        c = conn[j]
        for i in eachindex(c)
            ni = c[i]
            push!(map[ni], j)
        end
    end
    return map
end

"""
    FENodeToFEMap(conn::Vector{NTuple{N, IT}}, nmax::FInt) where {N, IT<:Integer}

Map from finite element nodes to the finite elements connecting them.

- `conns` = connectivities as a vector of tuples
- `nmax` = largest possible node number

Example:
```
m = FENodeToFEMap(fes.conn, count(fens))
```
"""
function FENodeToFEMap(conn::Vector{NTuple{N,IT}}, nmax::IT) where {N,IT<:Integer}
    return FENodeToFEMap(_makemap(conn, 1:length(conn), nmax))
end

"""
    FENodeToFEMap(fes::FE, nmax::IT) where {FE<:AbstractFESet,IT<:Integer}

Map from finite element nodes to the finite elements connecting them.

Convenience constructor.
"""
function FENodeToFEMap(fes::FE, nmax::IT) where {FE<:AbstractFESet,IT<:Integer}
    return FENodeToFEMap(_makemap(fes.conn, 1:count(fes), nmax))
end

"""
    FENodeToFEMap(fes::FE, nmax::IT) where {FE<:AbstractFESet,IT<:Integer}

Map from finite element nodes to the finite elements connecting them.

Convenience constructor.
"""
function FENodeToFEMap(mconn::Matrix{IT}, nmax::IT) where {IT<:Integer}
    map = Array{Vector{IT}}(undef, nmax)
    @inbounds for i in eachindex(map)
        map[i] = IT[]  # initially empty arrays
    end
    @inbounds for j in axes(mconn, 1)
        for i in axes(mconn, 2)
            push!(map[mconn[j, i]], j)
        end
    end
    return FENodeToFEMap(map)
end

end
