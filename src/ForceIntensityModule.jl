"""
    ForceIntensityModule

Module  to manage  distributed force intensity.
"""
module ForceIntensityModule

__precompile__(true)

import ..DataCacheModule: DataCache

"""
    ForceIntensity{T<:Number, F<:Function}

Distributed force (force intensity) type.

The force intensity class. The physical units are force per unit volume,
where volume depends on to which manifold the force is applied:
- force/length^3 (when applied to a 3-D solid),
- force/length^2 (when applied to a surface),
- force/length^1 (when applied along a curve), or
- force/length^0 (when applied at a point).

Signature of the function to compute the value of the force  at any given point
`XYZ`, using the columns of the Jacobian matrix of the element, `tangents`, the
finite element identifier, `feid`:

```
getforce!(forceout::Vector{CT}, XYZ::Matrix{T}, tangents::Matrix{T}, feid::IT) where {CT, T, IT}
```

A [`DataCache`](@ref) is used to store the data.
"""
struct ForceIntensity{DC<:DataCache}
    _cache::DC  # data cache  where the current value of the force can be retrieved
end

"""
    ForceIntensity(
        ::Type{T},
        ndofn,
        computeforce!::F,
    ) where {T<:Number, F<:Function}

Construct force intensity when the function to compute the intensity
vector is given.

# Arguments
- `T` = the type of the elements of the force vector, typically floating-point
  or complex floating-point numbers,
- `ndofn` = number of elements of the force vector (the length of the force
  vector),
- `computeforce!` = callback function. The function `computeforce!` needs to
  have a signature of
    ```
    function computeforce!(forceout::Vector{CT}, XYZ::Matrix{T},
        tangents::Matrix{T}, feid::IT) ) where {CT, T<:Number, IT<:Integer}
        # Calculate the force  and copy it into the buffer....
        return forceout
    end
    ```
    and it needs to  fill in the buffer `forceout` with the current force at the
    location `XYZ`, using, if appropriate, the information supplied in the Jacobian
    matrix `tangents`, the identifier of the finite element, `feid`.
"""
function ForceIntensity(
    ::Type{CT},
    ndofn,
    computeforce!::F,
) where {CT<:Number, F<:Function}
    # Allocate the buffer to be ready for the first call
    return ForceIntensity(DataCache(zeros(CT, ndofn), computeforce!))
end

"""
    ForceIntensity(force::Vector{T}) where {T<:Number}

Construct force intensity when the constant `force` vector is given.
"""
function ForceIntensity(force::Vector{CT}) where {CT<:Number}
    return ForceIntensity(DataCache(deepcopy(force)))
end

"""
    ForceIntensity(force::T) where {T<:Number}

Construct force intensity when the force is given as a scalar value.

The dimension of the force vector in this case is 1.
"""
function ForceIntensity(force::CT) where {CT<:Number}
    return ForceIntensity(CT[force])
end

"""
    updateforce!(self::ForceIntensity, XYZ::Matrix{T}, tangents::Matrix{T}, feid::IT, qpid::IT) where {T<:Number, IT<:Integer}

Update the force intensity vector.

Returns a vector (stored in the cache `self.cache`).
"""
function updateforce!(self::ForceIntensity, XYZ::Matrix{T}, tangents::Matrix{T}, feid::IT, qpid::IT) where {T<:Number, IT<:Integer}
    return self._cache(XYZ, tangents, feid, qpid)
end

end
