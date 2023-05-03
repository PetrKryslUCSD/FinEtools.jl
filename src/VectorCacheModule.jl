"""
    VectorCacheModule

Module to manage vector caches.
"""
module VectorCacheModule

__precompile__(true)

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
fillcache!(cacheout::Vector{T},
        XYZ::Matrix{T}, tangents::Matrix{T},
        fe_label; time::T = 0.0) where {T}
```
The cache `cacheout` is filled with the value  of the vector.
"""
struct VectorCache{CT<:Number, T<:Number, F<:Function}
    # Function to update and retrieve the vector
    _fillcache!::F
    # Cache where the current value of the vector can be retrieved
    _cache::Vector{CT}
    # Current time (or current load factor). Do not set directly. Use `settime!`
    _time::Ref{T}
end

"""
    VectorCache(::Type{CT}, nentries::IT, fillcache!::F) where {CT<:Number, F<:Function, IT}

Construct vector cache. The function to fill the vector cache is given.

This constructor is intended for *time-independent* vector caches.
This function needs to have a signature of
```
fillcache!(cacheout::Vector{CT},
        XYZ::Matrix{T}, tangents::Matrix{T},
        fe_label; time::T = 0.0) where {T}
    Calculate the vector and copy it into the cache....
    return forceout
end
```
and it needs to  fill in the cache `cacheout` with the current vector at
the location `XYZ`, using the information supplied in the Jacobian
matrix `tangents`, and the label of the finite element, `fe_label`, if
appropriate.
"""
function VectorCache(::Type{CT}, nentries::IT, fillcache!::F) where {CT<:Number, F<:Function, IT}
    T = real(CT(0.0))
    function fillcachenotime!(
        cacheout::Vector{CT},
        XYZ::Matrix{T},
        tangents::Matrix{T},
        fe_label;
        time::T = 0.0,
    ) where {CT, T}
        return fillcache!(cacheout, XYZ, tangents, fe_label)
    end
    # Allocate the cache to be ready for the first call
    return VectorCache(fillcachenotime!, zeros(CT, nentries), Ref(0.0))
end

"""
    VectorCache(
        ::Type{CT},
        nentries::IT,
        fillcache!::F,
        time::T,
    ) where {CT<:Number, F<:Function, IT}

Construct vector cache. The function to fill the vector cache is given.

This constructor is intended for *time-dependent* vector caches.
This function needs to have a signature of
```
fillcache!(cacheout::Vector{CT},
        XYZ::Matrix{T}, tangents::Matrix{T},
        fe_label; time::T = 0.0) where {CT, T}
    Calculate the vector and copy it into the cache....
    return forceout
end
```
and it needs to  fill in the cache `cacheout` with the current vector at
the location `XYZ`, using the information supplied in the Jacobian
matrix `tangents`, and the label of the finite element, `fe_label`, if
appropriate. The time can also be supplied (keyword argument `time`).
"""
function VectorCache(
    ::Type{CT},
    nentries::IT,
    fillcache!::F,
    time::T,
) where {CT<:Number, IT, F<:Function, T}
    # Allocate the cache to be ready for the first call
    return VectorCache(fillcache!, zeros(CT, nentries), Ref(time))
end

"""
    VectorCache(vector::Vector{CT}) where {CT<:Number}

Construct vector cache. The *constant* vector is given.
"""
function VectorCache(vector::Vector{CT}) where {CT<:Number}
    T = real(CT(0.0))
    function fillcache!(
        cacheout::Vector{CT},
        XYZ::Matrix{T},
        tangents::Matrix{T},
        fe_label;
        time::T = 0.0,
    ) where {CT, T}
        # do nothing:  the vector is already in the cache
        return cacheout
    end
    return VectorCache(fillcache!, deepcopy(vector), Ref(0.0))
end

"""
    updateretrieve!(self::VectorCache, XYZ::Matrix{T}, tangents::Matrix{T}, fe_label) where {T<:Number, IT}

Update the cache and retrieve the vector.

After the return from this function the updated vector can be read from
the cache as `self.cache` (also returned). If the vector depends on
time, the vector cache time first needs to be set as
```
settime!(c, t)
```
"""
function updateretrieve!(self::VectorCache, XYZ::Matrix{T}, tangents::Matrix{T}, fe_label) where {T<:Number}
    self._fillcache!(self._cache, XYZ, tangents, fe_label; time = self._time[])
    return self._cache
end

"""
    settime!(self::VectorCache, time)

Set the current time for the vector cache.
"""
function settime!(self::VectorCache, time)
    self._time[] = time
    return self
end

end
