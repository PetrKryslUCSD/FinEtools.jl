"""
    DataCacheModule

Module to manage array caches.
"""
module DataCacheModule

__precompile__(true)

import Base: size

"""
    DataCache{D, F<:Function, T<:Number}

Type for caching data, such as vectors, matrices, and numbers.

`D` = type of the data, for instance `Matrix{Float64}` or `Float32`.
`F` = type of the function to update the entries of the array.
`T` = type of the time,

Signature of the function to fill the cache with the value of the array is as
follows: it may use the location `XYZ`, it may use the columns of the Jacobian
matrix of the element, `tangents`, it may also choose the value given the
finite element label, `fe_label`. Finally, the value of the array may also
depend on the `time` (or the load factor). All of these values are chosen by
the code requesting the value of the cache. It must return the modified
argument `cacheout`.

```
function fillcache!(cacheout::D,
    XYZ::VecOrMat{T}, tangents::Matrix{T}, fe_label,
    time::T = 0.0) where {D, T}
    ... # modify the value of cacheout
    return cacheout
end
```

When the cache is accessed, the callback `fillcache!` is
called, and the output `cacheout` is filled with the value of the cached data.

# Example
```
function fillcache!(cacheout::Array{CT, N},
        XYZ::VecOrMat{T}, tangents::Matrix{T},
        fe_label, time::T = 0.0) where {CT, N, T}
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
data field is filled in and returned, but no data needs to be copied. The bad
news is, the cache is not thread safe. Reading is okay, but writing can lead to
data races.
"""
mutable struct DataCache{D, F<:Function}
    # Cache where the current value of the data can be retrieved
    _cache::D
    # Function to update and retrieve the array
    _fillcache!::F
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
        fe_label,
        time::T = 0.0,
    ) where {D, T}
        # do nothing:  the data is already in the cache
        return cacheout
    end
    return DataCache(deepcopy(data), _fillcache_constant!)
end

"""
    (c::DataCache)(XYZ::VecOrMat{T}, tangents::Matrix{T}, fe_label, time::T) where {T<:Number}

Update the cache and retrieve the array.
"""
function (c::DataCache)(XYZ::VecOrMat{T}, tangents::Matrix{T}, fe_label, time::T) where {T<:Number}
    c._cache = c._fillcache!(c._cache, XYZ, tangents, fe_label, time)
    return c._cache
end

"""
    (c::DataCache)(XYZ::VecOrMat{T}, tangents::Matrix{T}, fe_label) where {T<:Number}

Update the cache and retrieve the array.
"""
function (c::DataCache)(XYZ::VecOrMat{T}, tangents::Matrix{T}, fe_label) where {T<:Number}
    c._cache = c._fillcache!(c._cache, XYZ, tangents, fe_label, 0.0)
    return c._cache
end

"""
    size(self::DataCache)

Size of the data cache value.
"""
size(self::DataCache) = size(self._cache)

end # DataCacheModule
