"""
    SurfaceNormalModule

Module to evaluate surface normal vector.

The normal is assumed to be the exterior normal, and the vector is normalized to unit length.
"""
module SurfaceNormalModule

__precompile__(true)

import ..DataCacheModule: DataCache, updateretrieve!
using LinearAlgebra: cross, norm

"""
    SurfaceNormal{F<:Function}

Exterior surface normal type.

Normalized to unit length.

Signature of the function to compute the value of the unit normal
at any given point `XYZ`, using the columns of the Jacobian matrix
of the element, `tangents`, and if necessary  also the finite element label, `fe_label`:

```
computenormal!(normalout::Vector{T}, XYZ::Matrix{T}, tangents::Matrix{T}, fe_label<:Integer) where {T}
```

The buffer `normalout` is filled with the value  of the normal vector.
"""
struct SurfaceNormal{DC<:DataCache}
    # Cache of the current value of the normal
    _cache::DC
end

"""
    SurfaceNormal(ndimensions, computenormal!::F) where {F<:Function}

Construct surface normal evaluator when the function to compute the normal vector is
given. This function needs to have a signature of
```
function computenormal!(normalout::Vector{T}, XYZ::Matrix{T}, tangents::Matrix{T}, fe_label<:Integer) where {T}
    Calculate the normal and copy it into the buffer....
    return normalout # return the buffer
end
```
and it needs to  fill in the buffer `normalout` with the current vector at the
location `XYZ`, using if appropriate the information supplied in the Jacobian
matrix `tangents`, and the label of the finite element, `fe_label`.
"""
function SurfaceNormal(ndimensions, computenormal!::F) where {F<:Function}
    # Allocate the buffer to be ready for the first call
    return SurfaceNormal(DataCache(zeros(Float64, ndimensions), computenormal!))
end

"""
    SurfaceNormal(ndimensions, z::T, computenormal!::F) where {T<:Number, F<:Function}

Construct surface normal evaluator when the function to compute the normal vector is
given. This function needs to have a signature of
```
function computenormal!(normalout::Vector{T}, XYZ::Matrix{T}, tangents::Matrix{T}, fe_label<:Integer) where {T}
    Calculate the normal and copy it into the buffer....
    return normalout # return the buffer
end
```
and it needs to  fill in the buffer `normalout` with the current vector at the
location `XYZ`, using if appropriate the information supplied in the Jacobian
matrix `tangents`, and the label of the finite element, `fe_label`.
The type of the entries of the normal are `T`.
"""
function SurfaceNormal(ndimensions, z::T, computenormal!::F) where {T<:Number, F<:Function}
    # Allocate the buffer to be ready for the first call
    return SurfaceNormal(DataCache(zeros(T, ndimensions), computenormal!))
end

"""
    SurfaceNormal(ndimensions::FInt)

Construct surface normal evaluator when the default calculation of the normal
vector based on the columns of the Jacobian matrix should be used. 

The normal vector has `ndimensions` entries.

When the columns of the `tangents` array are parallel (or one of them is a zero
vector), the normal cannot be normalized to unit length (it is a zero vector).
In that case a zero vector is returned, and a warning is printed.
"""
function SurfaceNormal(ndimensions)
    function defaultcomputenormal!(
        normalout::Vector{T},
        XYZ::Matrix{T},
        tangents::Matrix{T},
        fe_label; time=0.0
    ) where {T}
        fill!(normalout, zero(T))
        # Produce a default normal
        if (size(tangents, 1) == 3) && (size(tangents, 2) == 2)# surface in three dimensions
            normalout[:] .= cross(vec(tangents[:, 1]), vec(tangents[:, 2]))# outer normal to the surface
        elseif (size(tangents, 1) == 2) && (size(tangents, 2) == 1)# curve in two dimensions
            normalout[1] = +tangents[2, 1]
            normalout[2] = -tangents[1, 1]# outer normal to the contour
        else
            error("No definition of normal vector")
        end
        nn = norm(normalout)
        if nn == 0.0
            @warn("Zero-length normal")
        else
            normalout ./= nn
        end
        return normalout
    end
    return SurfaceNormal(DataCache(zeros(Float64, ndimensions), defaultcomputenormal!))
end

"""
    SurfaceNormal(vector::Vector{T}) where {T<:Number}

Construct surface normal vector when the *constant* normal vector is given.
"""
function SurfaceNormal(vector::Vector{T}) where {T<:Number}
    return SurfaceNormal(DataCache(deepcopy(vector)))
end

"""
    updatenormal!(self::SurfaceNormal, XYZ::Matrix{T}, tangents::Matrix{T}, fe_label) where {T}

Update the surface normal vector.

Returns a vector (stored in the cache).
"""
function updatenormal!(self::SurfaceNormal, XYZ::Matrix{T}, tangents::Matrix{T}, fe_label) where {T}
    return self._cache(XYZ, tangents, fe_label)
end

end
