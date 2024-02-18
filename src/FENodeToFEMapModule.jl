"""
    FENodeToFEMapModule

Module to construct a map from finite element nodes to the finite elements.
"""
module FENodeToFEMapModule

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
    map = Vector{IT}[]
    sizehint!(map, nmax)
    for i = 1:nmax
        push!(map, [])  # initially empty arrays
    end
    for j in eachindex(conn)
        for i in eachindex(conn[j])
            ni = conn[j][i]
            push!(map[ni], j)
        end
    end
    return FENodeToFEMap(map)
end

"""
    FENodeToFEMap(conns::FIntMat, nmax::FInt)

Map from finite element nodes to the finite elements connecting them.

- `conns` = integer array of the connectivities
- `nmax` = largest possible node number
"""
function FENodeToFEMap(conns::Matrix{IT}, nmax::IT) where {IT<:Integer}
    map = Vector{IT}[]
    sizehint!(map, nmax)
    for i = 1:nmax
        push!(map, [])  # initially empty arrays
    end
    for j in axes(conns, 1)
        for i in axes(conns, 2)
            ni = conns[j, i]
            push!(map[ni], j)
        end
    end
    return FENodeToFEMap(map)
end

end
