"""
    ForceIntensityModule

Module  to manage  distributed force intensity.
"""
module ForceIntensityModule

__precompile__(true)

import ..DataCacheModule: DataCache, updateretrieve!, setcachetime!

"""
    ForceIntensity{T<:Number, F<:Function}

Distributed force (force intensity) type.

The force intensity class. The physical units are force per unit volume,
where volume depends on to which manifold the force is applied:
- force/length^3 (when applied to a 3-D solid),
- force/length^2 (when applied to a surface),
- force/length^1 (when applied along a curve), or
- force/length^0 (when applied at a point).

Signature of the function to compute the value of the force  at any
given point `XYZ`, using the columns of the Jacobian matrix of the
element, `tangents`, and if necessary  also the finite element label,
`fe_label`:
```
getforce!(forceout::Vector{CT}, XYZ::Matrix{T}, tangents::Matrix{T}, fe_label) where {CT, T}
```
The buffer `forceout` is filled with the value  of the force. The vector
`forceout` is also returned for convenience.
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

This constructor is intended for *time-independent* vector caches.

# Arguments
- `T` = the type of the elements of the force vector, typically floating-point or complex floating-point numbers,
- `ndofn` = number of elements of the force vector (the length of the force vector),
- `computeforce!` = callback function.
The function `computeforce!` needs to have a signature of
```
function computeforce!(forceout::Vector{CT}, XYZ::Matrix{T},
    tangents::Matrix{T}, fe_label) where {CT, T}
    # Calculate the force  and copy it into the buffer....
    return forceout
end
```
and it needs to  fill in the buffer `forceout` with the current force at the
location `XYZ`, using, if appropriate, the information supplied in the Jacobian
matrix `tangents`, and the label of the finite element, `fe_label`.
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
    ForceIntensity(
        ::Type{T},
        ndofn,
        computeforce!::F,
        time::T,
    ) where {T<:Number, F<:Function}

Construct force intensity when the function to compute the intensity
vector is given.

This constructor is intended for time-dependent force intensity caches.

# Arguments
- `T` = the type of the elements of the force vector, typically floating-point or complex floating-point numbers,
- `ndofn` = number of elements of the force vector (the length of the force vector),
- `computeforce!` = callback function,
- `time` = initial time.
The function `computeforce!` needs to have a signature of
```
function computeforce!(forceout::Vector{CT}, XYZ::Matrix{T},
    tangents::Matrix{T}, fe_label, time::T) where {CT, T}
	# Calculate the force  and copy it into the buffer....
	return forceout
end
```
and it needs to  fill in the buffer `forceout` with the current force at
the location `XYZ`, using if appropriate the information supplied in the
Jacobian matrix `tangents`, and the label of the finite element,
`fe_label`. The initial `time` is given.

The time needs to be set with `setcachetime!` before calling `updateforce!`
as follows:
```
XYZ = reshape([0.0, 0.0], 2, 1)
tangents = reshape([0.0, 1.0], 2, 1)
fe_label = 0
setvector!(v, XYZ, tangents, fe_label; time::FFlt = 0.0) = begin
    return (time < 5.0 ?  v .= [10.0] : v .= [0.0])
end
vector = [10.0]
fi = ForceIntensity(FFlt, length(vector), setvector!, 0.0)
v = updateforce!(fi, XYZ, tangents, fe_label)
@test v == [10.0]
setcachetime!(fi, 6.0)
v = updateforce!(fi, XYZ, tangents, fe_label)
@test v == [0.0]
```
"""
function ForceIntensity(
    ::Type{CT},
    ndofn,
    computeforce!::F,
    time::T,
) where {CT<:Number, F<:Function, T}
    # Allocate the buffer to be ready for the first call
    return ForceIntensity(setcachetime!(DataCache(zeros(CT, ndofn), computeforce!), time))
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
    updateforce!(self::ForceIntensity, XYZ::Matrix{T}, tangents::Matrix{T}, fe_label) where {T}

Update the force intensity vector.

Returns a vector (stored in the cache `self.cache`).
"""
function updateforce!(self::ForceIntensity, XYZ::Matrix{T}, tangents::Matrix{T}, fe_label) where {T}
    return self._cache(XYZ, tangents, fe_label)
end

"""
    setcachetime!(self::ForceIntensity, time)

Set the current time for the force intensity.
"""
function setcachetime!(self::ForceIntensity, time)
    setcachetime!(self._cache, time)
    return self
end

"""
    getcachetime(self::ForceIntensity)

Retrieve the current time for the force intensity.
"""
getcachetime(self::ForceIntensity) = getcachetime(self._cache)

end
