"""
    SurfaceNormalModule

Module to evaluate surface normal vector.

The normal is assumed to be the exterior normal, and the vector is normalized to unit length.
"""
module SurfaceNormalModule

__precompile__(true)

using ..FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import ..VectorCacheModule: VectorCache, updateretrieve!
using LinearAlgebra: cross, norm

"""
    SurfaceNormal{F<:Function}

Exterior surface normal type.

Normalized to unit length.

Signature of the function to compute the value of the unit normal
at any given point `XYZ`, using the columns of the Jacobian matrix
of the element, `tangents`, and if necessary  also the finite element label, `fe_label`:

```
computenormal!(normalout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
```

The buffer `normalout` is filled with the value  of the normal vector.
"""
struct SurfaceNormal{F<:Function}
	"""
	    Cache of the current value of the normal
	"""
    cache::VectorCache{FFlt, F}
end

"""
    SurfaceNormal(ndimensions::FInt, computenormal!::F) where {F<:Function}

Construct surface normal evaluator when the function to compute the normal vector is
given. This function needs to have a signature of
```
function computenormal!(normalout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    Calculate the normal and copy it into the buffer....
    return normalout # return the buffer
end
```
and it needs to  fill in the buffer `normalout` with the current vector at the
location `XYZ`, using if appropriate the information supplied in the Jacobian
matrix `tangents`, and the label of the finite element, `fe_label`.
"""
function SurfaceNormal(ndimensions::FInt, computenormal!::F) where {F<:Function}
    # Allocate the buffer to be ready for the first call
    return SurfaceNormal(VectorCache(FFlt, ndimensions, computenormal!));
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
function SurfaceNormal(ndimensions::FInt)
    function defaultcomputenormal!(normalout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        @assert !any(x->isnan(x), tangents[:])
        fill!(normalout, 0.0)
        # Produce a default normal
        if (size(tangents,1) == 3) && (size(tangents,2) == 2)# surface in three dimensions
            normalout[:] .= cross(vec(tangents[:,1]),vec(tangents[:,2]));# outer normal to the surface
        elseif (size(tangents,1)==2)  && (size(tangents,2)==1)# curve in two dimensions
            normalout[1] = +tangents[2,1];
            normalout[2] = -tangents[1,1];# outer normal to the contour
        else
            error("No definition of normal vector");
        end
        @assert !any(x->isnan(x), normalout)
        nn = norm(normalout);
        if  nn == 0.0 
            @warn("Zero-length normal")
        else
            normalout ./= nn
        end
        return normalout
    end
    return SurfaceNormal(VectorCache(FFlt, ndimensions, defaultcomputenormal!));
end

"""
    SurfaceNormal(vector::FVec{T}) where {T<:Number}

Construct surface normal vector when the *constant* normal vector is given.
"""
function SurfaceNormal(vector::FVec{T}) where {T<:Number}
    return SurfaceNormal(VectorCache(deepcopy(vector)));
end

"""
    updatenormal!(self::SurfaceNormal, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)

Update the surface normal vector.

Returns a vector (stored in the cache).
"""
function updatenormal!(self::SurfaceNormal, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    return updateretrieve!(self.cache, XYZ, tangents, fe_label)
end

end
