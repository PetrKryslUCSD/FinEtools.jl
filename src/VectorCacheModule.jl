"""
    VectorCacheModule

Module to manage vector caches.
"""
module VectorCacheModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict

"""
    VectorCache{T<:Number, F<:Function}

Type for caching vectors.

`T` = type of the entries of the vector,
`F` = type of the function to update the entries of the vector.

Signature of the function to fill the cache with the value of the vector at any given point
`XYZ`, using the columns of the Jacobian matrix of the element, `tangents`,
and if necessary  also the finite element label, `fe_label`:
```
fillcache!(cacheout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
```
The cache `cacheout` is filled with the value  of the vector.
"""
struct VectorCache{T<:Number, F<:Function}
    fillcache!::F # Function to update and retrieve the vector
    cache::Vector{T}    # Cache where the current value of the vector can be retrieved
end

"""
    VectorCache(::Type{T}, nentries::FInt, fillcache!::F) where {T<:Number, F<:Function}

Construct vector cache. The function to fill the vector cache is
given. This function needs to have a signature of
```
fillcache!(cacheout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
Calculate the vector and copy it into the cache....
return forceout
end
```
and it needs to  fill in the cache `cacheout` with the current vector at the
location `XYZ`, using the information supplied in the Jacobian
matrix `tangents`, and the label of the finite element, `fe_label`, if appropriate.
"""
function VectorCache(::Type{T}, nentries::FInt, fillcache!::F) where {T<:Number, F<:Function}
    # Allocate the cache to be ready for the first call
    return VectorCache(fillcache!, zeros(T, nentries));
end

"""
    VectorCache(vector::FVec{T}) where {T<:Number}

Construct vector cache. The constant vector is given.
"""
function VectorCache(vector::FVec{T}) where {T<:Number}
    function fillcache!(cacheout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        # do nothing:  the vector is already in the cache
        return cacheout
    end
    return VectorCache(fillcache!, deepcopy(vector));
end

"""
    updateretrieve!(self::VectorCache, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)

Update the cache and retrieve the vector.

After the return from this function the updated vector can be read from the
cache as `self.cache` (also returned).
"""
function updateretrieve!(self::VectorCache, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    self.fillcache!(self.cache, XYZ, tangents, fe_label)
    return self.cache
end

end
