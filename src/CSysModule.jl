"""
    CSysModule

Module for management of coordinate systems.
"""
module CSysModule

__precompile__(true)

import LinearAlgebra: norm, cross
using ..RotationUtilModule: cross3!

"""
    CSys{T<:Number, F<:Function}

Type for coordinate system transformations. Used to define material coordinate
systems, and output coordinate systems, for instance.
"""
struct CSys{T <: Number, F <: Function}
    isconstant::Bool
    isidentity::Bool
    __updatebuffer!::F # function to update the coordinate system matrix.
    # Signature: `update!(csmatout::Matrix{T}, XYZ::VecOrMat{T}, tangents::Matrix{T},
    # feid::IT, qpid::IT) where {T, IT}`
    _csmat::Array{T, 2} # the coordinate system matrix (buffer); see
end

"""
    CSys(sdim::IT1, mdim::IT2, computecsmat::F) where {F <: Function, IT1, IT2}

Construct coordinate system when the function to compute the
rotation matrix is given.

The function signature:
```
update!(csmatout::Matrix{T}, XYZ::VecOrMat{T}, tangents::Matrix{T},
    feid::IT, qpid::IT) where {T, IT}
```
where
- `csmatout`= output matrix buffer, of size `(sdim, mdim)`
- `XYZ`= location  in physical coordinates,
- `tangents`= tangent vector matrix, tangents to the parametric coordinate
  curves  in the element,
- `feid`= finite element identifier;
- `qpid`= quadrature point identifier.

# Example

```
# Cylindrical coordinate system: NO ALLOCATIONS WHATSOEVER!
@views function compute!(csmatout, XYZ, tangents, feid, qpid)
    center = (0.0, 0.0, 0.0)
    xyz = (XYZ[1], XYZ[2], XYZ[3])
    csmatout[:, 1] .= xyz .- center
    csmatout[3, 1] = 0.0
    csmatout[:, 1] ./= norm(csmatout[:, 1])
    csmatout[:, 3] .= (0.0, 0.0, 1.0)
    cross3!(csmatout[:, 2], csmatout[:, 3], csmatout[:, 1])
    csmatout[:, 2] ./=  norm(csmatout[:, 2])
    return csmatout
end
```
"""
function CSys(sdim::IT1, mdim::IT2, computecsmat::F) where {F <: Function, IT1, IT2}
    csmat = fill(zero(Float64), sdim, mdim) # Allocate buffer for the first call
    return CSys(false, false, computecsmat, csmat)
end

"""
    CSys(sdim::IT1, mdim::IT2, z::T, computecsmat::F) where {IT1, IT2, T <: Number, F <: Function}

Construct coordinate system when the function to compute the
rotation matrix of type `T` is given.


- `z` = zero value,
- The `computecsmat` function signature:
    ```
    update!(csmatout::Matrix{T}, XYZ::VecOrMat{T}, tangents::Matrix{T},
        feid::IT, qpid::IT) where {T, IT}
    ```
    where
    - `csmatout`= output matrix buffer, of size `(sdim, mdim)`;
    - `XYZ`= location  in physical coordinates;
    - `tangents`= tangent vector matrix, tangents to the parametric coordinate
      curves  in the element;
    - `feid`= finite element identifier;
    - `qpid`= quadrature point identifier.


# Example

```
# Cylindrical coordinate system: NO ALLOCATIONS WHATSOEVER!
@views function compute!(csmatout, XYZ, tangents, feid, qpid)
    center = (0.0, 0.0, 0.0)
    xyz = (XYZ[1], XYZ[2], XYZ[3])
    csmatout[:, 1] .= xyz .- center
    csmatout[3, 1] = 0.0
    csmatout[:, 1] ./= norm(csmatout[:, 1])
    csmatout[:, 3] .= (0.0, 0.0, 1.0)
    cross3!(csmatout[:, 2], csmatout[:, 3], csmatout[:, 1])
    csmatout[:, 2] ./=  norm(csmatout[:, 2])
    return csmatout
end
```
"""
function CSys(sdim::IT1, mdim::IT2, z::T, computecsmat::F) where {IT1, IT2, T <: Number, F <: Function}
    csmat = fill(z, sdim, mdim) # Allocate buffer, in preparation for the first call
    return CSys(false, false, computecsmat, csmat)
end

"""
    CSys(csmat::Matrix{T}) where {T}

Construct coordinate system when the rotation matrix is given.
"""
function CSys(csmat::Matrix{T}) where {T}
    function __updatebuffer!(csmatout::Matrix{T},
        XYZ::VecOrMat{T},
        tangents::Matrix{T},
        feid::IT1, qpid::IT2) where {T, IT1, IT2}
        return csmatout # nothing to be done here, the matrix is already in the buffer
    end
    return CSys(true, false, __updatebuffer!, deepcopy(csmat))# fill the buffer with the given matrix
end

"""
    CSys(dim, z::T) where {T}

Construct coordinate system when the rotation matrix of element type `T` is the
identity.

`dim` = is the space dimension.
"""
function CSys(dim::IT, z::T) where {IT<:Integer, T}
    function __updatebuffer!(csmatout::Matrix{T},
        XYZ::VecOrMat{T},
        tangents::Matrix{T},
        feid::IT1, qpid::IT2) where {T, IT1, IT2}
        return csmatout # nothing to be done here, the matrix is already in the buffer
    end
    return CSys(true,
        true,
        __updatebuffer!,
        [i == j ? one(T) : zero(T) for i in 1:dim, j in 1:dim])# identity
end

"""
    CSys(dim::IT) where {IT}

Construct coordinate system when the rotation matrix is the identity.

`dim` = is the space dimension.
"""
function CSys(dim::IT) where {IT}
    return CSys(dim, Float64(0.0))
end

"""
    CSys(sdim::IT1, mdim::IT2) where {IT1<:Integer, IT2<:Integer}

Construct coordinate system for isotropic-material used with isoparametric
finite elements.

- `sdim` = number of space dimensions,
- `mdim` = number of manifold dimensions of the finite element in which the
  coordinate system  is being evaluated.

!!! note

    If  the coordinate system matrix should be identity, better use the
    constructor for this specific situation, `CSys(dim)`. That will be much
    more efficient.

# See also
`gen_iso_csmat`
"""
function CSys(sdim::IT1, mdim::IT2) where {IT1<:Integer, IT2<:Integer}
    function __updatebuffer!(csmatout::Matrix{T},
        XYZ::VecOrMat{T},
        tangents::Matrix{T},
        feid::IT1, qpid::IT2) where {T, IT1, IT2}
        return gen_iso_csmat!(csmatout, XYZ, tangents, feid, qpid)
    end
    return CSys(false, false, __updatebuffer!, fill(zero(Float64), sdim, mdim))
end

"""
    updatecsmat!(self::CSys,
        XYZ::Matrix{T},
        tangents::Matrix{T},
        feid::IT1,
        qpid::IT2) where {T, IT1, IT2}

Update the coordinate system orientation matrix.

The  coordinate system matrix is updated based upon the location `XYZ` of the
evaluation point, and possibly on the Jacobian matrix `tangents` within the
element in which the coordinate system matrix is evaluated,  or perhaps on the
identifier `feid` of the finite element and/or the quadrature point identifier.

After this function returns, the coordinate system matrix can be read in the
buffer as `self.csmat`.
"""
function updatecsmat!(self::CSys,
    XYZ::Matrix{T},
    tangents::Matrix{T},
    feid::IT1,
    qpid::IT2) where {T, IT1, IT2}
    self.__updatebuffer!(self._csmat, XYZ, tangents, feid, qpid)
    return self._csmat
end

"""
    csmat(self::CSys)

Return coordinate system rotation matrix.

No allocation is involved.
"""
function csmat(self::CSys)
    self._csmat
end

"""
    gen_iso_csmat!(csmatout::Matrix{T}, XYZ::Matrix{T}, tangents::Matrix{T}, feid::IT, qpid::IT) where {T, IT}

Compute the coordinate system  for an isotropic material using information
available  by looking at the coordinate curves of isoparametric finite
elements.

- `XYZ`= location  in physical coordinates,
- `tangents`= tangent vector matrix, tangents to the parametric coordinate
  curves  in the element,
- `feid`= finite element identifier;
- `qpid` = quadrature point identifier.

The basic assumption here is that the material is isotropic, and
therefore the choice of the material directions does not really matter as
long as they correspond to the dimensionality of the element. For
instance a one-dimensional element (L2 as an example) may be embedded
in a three-dimensional space.

This function assumes that it is being called for an `mdim`-dimensional manifold
element, which is embedded in a `sdim`-dimensional Euclidean space. If `mdim ==
sdim`, the coordinate system matrix is the identity; otherwise the local
coordinate directions are aligned with the linear subspace defined by the
tangent vectors.

!!! warning

    This *cannot* be reliably used to produce consistent stresses because each
    quadrature point gets a local coordinate system which depends on the
    orientation of the element, in general different from the neighboring elements.
"""
@views function gen_iso_csmat!(csmatout::Matrix{T},
    XYZ::Matrix{T},
    tangents::Matrix{T},
    feid::IT1,
    qpid::IT2) where {T, IT1, IT2}
    sdim, mdim = size(tangents)
    if sdim == mdim # finite element embedded in space of the same dimension
        for i in 1:size(csmatout, 1), j in 1:size(csmatout, 2)
            csmatout[i, j] = (i == j ? one(T) : zero(T))
        end
    else # lower-dimensional finite element embedded in space of higher dimension
        @assert 0 < mdim < sdim
        @assert 0 < sdim
        csmatout[:, 1] = tangents[:, 1]
        csmatout[:, 1] ./= norm(csmatout[:, 1])
        if mdim == 1 # curve-like finite element in 2d or 3d
        # all done
        elseif mdim == 2 # surface-like finite element in 3d
            e2 = (tangents[1, 2], tangents[2, 2], tangents[3, 2])
            cross3!(csmatout[:, 2], csmatout[:, 1], e2)
            e3 = (csmatout[1, 2], csmatout[2, 2], csmatout[3, 2])
            cross3!(csmatout[:, 2], e3, csmatout[:, 1])
            csmatout[:, 2] ./= norm(csmatout[:, 2])
        end
    end
    return csmatout
end

end # module
