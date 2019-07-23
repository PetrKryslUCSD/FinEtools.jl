"""
    VectorCacheModule

Module to manage vector caches.
"""
module VectorCacheModule

using ..FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict

"""
    VectorCache{T<:Number, F<:Function}

Type for caching vectors.

`T` = type of the entries of the vector,
`F` = type of the function to update the entries of the vector.

Signature of the function to fill the cache with the value of the vector
at any given point `XYZ`, using the columns of the Jacobian matrix of
the element, `tangents`, and, if convenient, also the finite element
label, `fe_label`. Finally, the value of the vector may also depend on
the `time` (or the load factor):
```
fillcache!(cacheout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt; time::FFlt = 0.0)
```
The cache `cacheout` is filled with the value  of the vector.
"""
struct VectorCache{T<:Number, F<:Function}
	# Function to update and retrieve the vector
	fillcache!::F
    # Cache where the current value of the vector can be retrieved
    cache::Vector{T}
	# Current time (or current load factor). Do not set directly. Use `settime!`
	_time::Vector{FFlt}
end

"""
    VectorCache(::Type{T}, nentries::FInt, fillcache!::F) where {T<:Number, F<:Function}

Construct vector cache. The function to fill the vector cache is given.

This constructor is intended for *time-independent* vector caches.
This function needs to have a signature of
```
fillcache!(cacheout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    Calculate the vector and copy it into the cache....
    return forceout
end
```
and it needs to  fill in the cache `cacheout` with the current vector at
the location `XYZ`, using the information supplied in the Jacobian
matrix `tangents`, and the label of the finite element, `fe_label`, if
appropriate.
"""
function VectorCache(::Type{T}, nentries::FInt, fillcache!::F) where {T<:Number, F<:Function}
	function fillcachenotime!(cacheout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt; time::FFlt = 0.0)
        return fillcache!(cacheout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    end
    # Allocate the cache to be ready for the first call
    return VectorCache(fillcachenotime!, zeros(T, nentries), [0.0]);
end

"""
    VectorCache(::Type{T}, nentries::FInt, fillcache!::F, time::FFlt) where {T<:Number, F<:Function}

Construct vector cache. The function to fill the vector cache is given.

This constructor is intended for *time-dependent* vector caches.
This function needs to have a signature of
```
fillcache!(cacheout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt; time::FFlt = 0.0)
    Calculate the vector and copy it into the cache....
    return forceout
end
```
and it needs to  fill in the cache `cacheout` with the current vector at
the location `XYZ`, using the information supplied in the Jacobian
matrix `tangents`, and the label of the finite element, `fe_label`, if
appropriate. The time can also be supplied (keyword argument `time`).
"""
function VectorCache(::Type{T}, nentries::FInt, fillcache!::F, time::FFlt) where {T<:Number, F<:Function}
	# Allocate the cache to be ready for the first call
    return VectorCache(fillcache!, zeros(T, nentries), [time]);
end

"""
    VectorCache(vector::FVec{T}) where {T<:Number}

Construct vector cache. The *constant* vector is given.
"""
function VectorCache(vector::FVec{T}) where {T<:Number}
    function fillcache!(cacheout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt; time::FFlt = 0.0)
        # do nothing:  the vector is already in the cache
        return cacheout
    end
    return VectorCache(fillcache!, deepcopy(vector), [0.0]);
end

"""
    updateretrieve!(self::VectorCache, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)

Update the cache and retrieve the vector.

After the return from this function the updated vector can be read from
the cache as `self.cache` (also returned). If the vector depends on
time, the vector cache time first needs to be set as
```
settime!(c, t)
```
"""
function updateretrieve!(self::VectorCache, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    self.fillcache!(self.cache, XYZ, tangents, fe_label; time = self._time[1])
    return self.cache
end

"""
    settime!(self::VectorCache, time::FFlt)

Set the current time for the vector cache.
"""
function settime!(self::VectorCache, time::FFlt)
	self._time[1] = time
	return self
end

end
