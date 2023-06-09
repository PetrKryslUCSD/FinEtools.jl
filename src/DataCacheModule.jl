"""
    DataCacheModule

Module to manage array caches.
"""
module DataCacheModule

__precompile__(true)

"""
    DataCache{D, F<:Function, T<:Number}

Type for caching data, such as vectors, matrices, and numbers.

`D` = type of the data, for instance `Matrix{Float64}` or `Float32`.
`F` = type of the function to update the entries of the array.
`T` = type of the time,

Signature of the function to fill the cache with the value of the array,
at any given point `XYZ`, using the columns of the Jacobian matrix of
the element, `tangents`, and, if convenient, also the finite element
label, `fe_label`. Finally, the value of the array may also depend on
the `time` (or the load factor):
```
function fillcache!(cacheout::D,
        XYZ::VecOrMat{T}, tangents::Matrix{T}, fe_label;
        time::T = 0.0) where {D, T}
        ...
end
```!
When the cache is accessed with `updateretrieve!`, the callback `fillcache!` is
called, and the output `cacheout` is filled with the value of the cached data.

# Example
```
function fillcache!(cacheout::Array{CT, N},
        XYZ::VecOrMat{T}, tangents::Matrix{T},
        fe_label; time::T = 0.0) where {CT, N, T}
    cacheout .= LinearAlgebra.I(3)
    return cacheout
end
c = DataCache(zeros(Float32, 3, 3), fillcache!)
function f(c)
    XYZ, tangents, fe_label = (reshape([0.0, 0.0], 1, 2), [1.0 0.0; 0.0 1.0], 1)
    data = c(XYZ, tangents, fe_label)
end
@test f(c) == LinearAlgebra.I(3)
```

!!! note

The point of the data cache is that there will be no copying of data. The cache
data field is filled in and returned, but no data needs to be copied. The bad news
is, the cache is not thread safe. Reading is okay, but writing can lead to data
races.
"""
mutable struct DataCache{D, F<:Function, T<:Number}
    # Function to update and retrieve the array
    _fillcache!::F
    # Cache where the current value of the data can be retrieved
    _cache::D
    # Current time (or current load factor). Do not set directly. Use `setcachetime!`
    _time::Ref{T}
end

"""
    DataCache(data::D, fillcache!::F) where {D, F<:Function}

Construct data cache. The function to fill the data in the cache is given.

The callback function needs to have a signature of
```
function fillcache!(cacheout::D,
        XYZ::VecOrMat{T}, tangents::Matrix{T},
        fe_label; time::T = 0.0) where {D, T}
    #   Calculate the data and copy it into the cacheout data
    # or simply return it if it is a number value....
    return cacheout
end
```
and it needs to  fill in the cache argument `cacheout` with the current value at
the location `XYZ`, using the information supplied in the Jacobian
matrix `tangents`, and the label of the finite element, `fe_label`, if
appropriate. Time can also be taken into account.
"""
function DataCache(data::D, fillcache!::F) where {D, F<:Function}
    function _fillcache!(
        cacheout::D,
        XYZ::VecOrMat{T},
        tangents::Matrix{T},
        fe_label;
        time::T = 0.0
    ) where {D, T}
        return fillcache!(cacheout, XYZ, tangents, fe_label; time=time)
    end
    return DataCache(_fillcache!, deepcopy(data), Ref(0.0))
end

"""
    DataCache(data::Array{CT, N}) where {CT<:Number, N}

Construct data cache. The *constant* data is given.
"""
function DataCache(data::D) where {D}
    function _fillcache_constant!(
        cacheout::D,
        XYZ::VecOrMat{T},
        tangents::Matrix{T},
        fe_label;
        time::T = 0.0,
    ) where {D, T}
        # do nothing:  the data is already in the cache
        return cacheout
    end
    return DataCache(_fillcache_constant!, deepcopy(data), Ref(0.0))
end

"""
    updateretrieve!(self::C, XYZ::VecOrMat{T}, tangents::Matrix{T}, fe_label) where {C<:DataCache, T<:Number}

Update the cache and retrieve the array.

This function returns the cached array, filled with the values provided by the
callback.

If the array depends on time, the array cache time first needs to be set as
```
setcachetime!(c, t)
```
"""
function updateretrieve!(self::C, XYZ::VecOrMat{T}, tangents::Matrix{T}, fe_label) where {C<:DataCache, T<:Number}
    self._cache = self._fillcache!(self._cache, XYZ, tangents, fe_label; time = self._time[])
    return self._cache
end

"""
    (c::DataCache)(XYZ::VecOrMat{T}, tangents::Matrix{T}, fe_label) where {T<:Number}

Update the cache and retrieve the array.

This function returns the cached array, filled with the values provided by the
callback.

If the array depends on time, the array cache time first needs to be set as
```
setcachetime!(c, t)
```
"""
function (c::DataCache)(XYZ::VecOrMat{T}, tangents::Matrix{T}, fe_label) where {T<:Number}
    c._cache = c._fillcache!(c._cache, XYZ, tangents, fe_label; time = c._time[])
    return c._cache
end

"""
    setcachetime!(self::DataCache, time)

Set the current time for the array cache.
"""
function setcachetime!(self::DataCache, currenttime)
    self._time[] = currenttime
    return self
end

"""
    getcachetime(self::DataCache)

Retrieved the current time for the array cache.
"""
getcachetime(self::DataCache) = self._time[]

end # DataCacheModule
