"""
    MeshUtilModule

Module for mesh utility functions used in other meshing modules.
"""
module MeshUtilModule

using ..FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import Base.BitSet
import LinearAlgebra: norm

"""
    HyperFaceContainer

Type of a hyper face data to store a newly generated node.
"""
mutable struct HyperFaceContainer
    o::BitSet # numbers of the other nodes on the hyperface
    n::Int # new node number generated on the hyperface
end

"""
    makecontainer()

Make hyper face container.

This is a dictionary of hyper faces, indexed with an anchor node. The anchor
node is the smallest node number within the connectivity of the hyper face.
"""
makecontainer() = Dict{FInt, Array{HyperFaceContainer}}();

"""
    addhyperface!(container,hyperface,newn)

Add a hyper face to the container.

The new node is stored in hyper face data in the container and can be retrieved
later.
"""
function addhyperface!(container,hyperface,newn)
    h=sort([i for i in hyperface])
    anchor=h[1]; other= BitSet(h[2:end]);
    C=get(container,anchor,HyperFaceContainer[]);
    fnd=false;
    for k=1:length(C)
        if C[k].o == other
            fnd=true; break;
        end
    end
    if fnd!=true
        push!(C,HyperFaceContainer(other,newn)); newn=newn+1;#new node added
    end
    container[anchor]= C; # set the new array
    return newn # return the number of the added node, possibly unchanged
end

"""
    findhyperface!(container,hyperface)

Find a hyper face in the container.

If the container holds the indicated hyper face, it is returned together with
the stored new node number.
"""
function findhyperface!(container,hyperface)
    h=sort([i for i in hyperface])
    anchor=h[1]; other= BitSet(h[2:end]);
    C=get(container,anchor,HyperFaceContainer[]);
    for k=1:length(C)
        if C[k].o == other
            return [anchor, C[k].o], C[k].n
        end
    end
    return [], 0
end

function ontosphere(xyz::FFltMat,radius::FFlt)
# % Project nodes onto a sphere of given radius.
# %
# % function fens= onto_sphere(fens,radius,list)
# %
# % fens = finite element node set,
# % radius = radius of the sphere,
# % list = optional argument, if not empty then only
# %	     the nodes in the list are to be moved; otherwise all nodes are moved.

    for j=1:size(xyz,1)
        xyz[j,:] =xyz[j,:]*radius/norm(xyz[j,:]);
    end
    return xyz;
end

"""
    linearspace(start::T, stop::T, N::Int)  where {T<:Number}

Generate linear space.

Generate a linear sequence of `N` numbers between  `start` and `stop` (i. e.
sequence  of number with uniform intervals inbetween).

# Example
```
julia> linearspace(2.0, 3.0, 5)
2.0:0.25:3.0
```
"""
function linearspace(start::T, stop::T, N::Int)  where {T<:Number}
    return range(start, stop = stop, length = N)
end

"""
    gradedspace(start::T, stop::T, N::Int)  where {T<:Number}

Generate quadratic space.

Generate a quadratic sequence of `N` numbers between `start` and `stop`. This
sequence corresponds to separation of adjacent numbers that increases linearly
from start to finish.

# Example
```
julia> gradedspace(2.0, 3.0, 5)
5-element Array{Float64,1}:
 2.0
 2.0625
 2.25
 2.5625
 3.0
```
"""
function gradedspace(start::T, stop::T, N::Int, strength=2)  where {T<:Number}
    x = range(0.0, stop = 1.0, length = N);
    x = x.^strength
    # for i = 1:strength
    #     x = cumsum(x);
    # end
    x = x/maximum(x);
    out = start .* (1.0 .- x) .+ stop .* x;
end

# function inbox(box::FFltVec,sdim::FInt,x::FFltVec)
#     for i=1:sdim
#         if (!inrange(box[2*i-1],box[2*i],x[i]))
#             return false
#         end
#     end
#     return true
# end

end
