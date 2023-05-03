"""
    CSysModule

Module for management of coordinate systems.
"""
module CSysModule

__precompile__(true)

import LinearAlgebra: norm, cross

"""
    CSys{T<:Number, F<:Function}

Type for coordinate system transformations. Used to define material coordinate
systems, and output coordinate systems, for instance.
"""
struct CSys{T<:Number, F<:Function}
    isconstant::Bool
    isidentity::Bool
    updatebuffer!::F # function to update the coordinate system matrix.
    # `updatebuffer!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)`
    _csmat::Array{T,2} # the coordinate system matrix (buffer); see
end


"""
    CSys(sdim, mdim, computecsmat::F) where {F<:Function}

Construct coordinate system when the function to compute the
rotation matrix is given.

The function signature:
```
update!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
```
where
- `csmatout`= output matrix buffer, of size `(sdim, mdim)`
- `XYZ`= location  in physical coordinates,
- `tangents`= tangent vector matrix, tangents to the parametric coordinate
  curves  in the element,
- `fe_label`= finite element label.
"""
function CSys(sdim, mdim, computecsmat::F) where {F<:Function}
    csmat = fill(zero(Float64), sdim, mdim) # Allocate buffer, in preparation for the first call
    return CSys(false, false, computecsmat, csmat)
end

"""
    CSys(sdim, mdim, z::T, computecsmat::F) where {T<:Number, F<:Function}

Construct coordinate system when the function to compute the
rotation matrix of type `T` is given.


- `z` = zero value,
- The `computecsmat` function signature:
    ```
    update!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    ```
    where
    - `csmatout`= output matrix buffer, of size `(sdim, mdim)`
    - `XYZ`= location  in physical coordinates,
    - `tangents`= tangent vector matrix, tangents to the parametric coordinate
      curves  in the element,
    - `fe_label`= finite element label.
"""
function CSys(sdim, mdim, z::T, computecsmat::F) where {T<:Number, F<:Function}
    csmat = fill(z, sdim, mdim) # Allocate buffer, in preparation for the first call
    return CSys(false, false, computecsmat, csmat)
end

"""
    CSys(csmat::Matrix{T}) where {T}

Construct coordinate system when the rotation matrix is given.
"""
function CSys(csmat::Matrix{T}) where {T}
    function updatebuffer!(
        csmatout::Matrix{T},
        XYZ::Matrix{T},
        tangents::Matrix{T},
        fe_label,
    )
        return csmatout # nothing to be done here, the matrix is already in the buffer
    end
    return CSys(true, false, updatebuffer!, deepcopy(csmat))# fill the buffer with the given matrix
end

"""
    CSys(dim, z::T) where {T}

Construct coordinate system when the rotation matrix of element type `T` is the
identity.

`dim` = is the space dimension.
"""
function CSys(dim, z::T) where {T}
    return CSys([i == j ? one(T) : zero(T) for i in 1:dim, j in 1:dim])
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
    CSys(sdim::IT, mdim::IT) where {IT}

Construct coordinate system for isotropic-material used with isoparametric
finite elements.

- `sdim` = number of space dimensions,
- `mdim` = number of manifold dimensions of the finite element in which the
  coordinate system  is being evaluated.

!!! note
    If  the coordinate system matrix  should be identity, better use the constructor
    for this specific situation, `CSys(dim)`. That will be much more efficient.

# See also
`gen_iso_csmat`
"""
function CSys(sdim::IT, mdim::IT) where {IT}
    csmat = fill(zero(Float64), sdim, mdim) # Allocate buffer, prepare for the first call
    function updatebuffer!(
        csmatout::Matrix{T},
        XYZ::Matrix{T},
        tangents::Matrix{T},
        fe_label,
    ) where {T}
        gen_iso_csmat!(csmatout, XYZ, tangents, fe_label)
        return csmatout
    end
    return CSys(false, false, updatebuffer!, csmat)
end

"""
    updatecsmat!(self::CSys, XYZ::Matrix{T}, tangents::Matrix{T}, fe_label) where {T}

Update the coordinate system orientation matrix.

The  coordinate system matrix is updated based upon the location `XYZ` of the
evaluation point, and possibly on the Jacobian matrix `tangents` within the
element in which the coordinate system matrix is evaluated,  or perhaps on the
label `fe_label` of the finite element.

After this function returns, the coordinate system matrix can be retrieved
from the buffer as `self.csmat`.
"""
function updatecsmat!(self::CSys, XYZ::Matrix{T}, tangents::Matrix{T}, fe_label) where {T}
    self.updatebuffer!(self._csmat, XYZ, tangents, fe_label)
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
    gen_iso_csmat!(csmatout::Matrix{T}, XYZ::Matrix{T}, tangents::Matrix{T}, fe_label) where {T}

Compute the coordinate system  for an isotropic material using information
available  by looking at the coordinate curves of isoparametric finite
elements.

- `XYZ`= location  in physical coordinates,
- `tangents`= tangent vector matrix, tangents to the parametric coordinate
  curves  in the element,
- `fe_label`= finite element label.

The basic assumption here is that the material is isotropic, and
therefore the choice of the material directions does not really matter as
long as they correspond to the dimensionality of the element. For
instance a one-dimensional element (L2 as an example) may be embedded
in a three-dimensional space.

This function assumes that it is being called for an mdim-dimensional manifold
element, which is embedded in a sdim-dimensional Euclidean space. If `mdim ==
sdim`, the coordinate system matrix is the identity; otherwise the local
coordinate directions are aligned with the linear subspace defined by the
tangent vectors.

!!! warning

This *cannot* be reliably used to produce consistent stresses because each
quadrature point gets a local coordinate system which depends on the
orientation of the element.
"""
function gen_iso_csmat!(csmatout::Matrix{T}, XYZ::Matrix{T}, tangents::Matrix{T}, fe_label) where {T}
    sdim, mdim = size(tangents)
    if sdim == mdim # finite element embedded in space of the same dimension
        copyto!(csmatout, [i == j ? one(T) : zero(T) for i in 1:sdim, j in 1:sdim])
    else # lower-dimensional finite element embedded in space of higher dimension
        @assert 0 < mdim < 3
        e1 = tangents[:, 1] / norm(tangents[:, 1])
        if mdim == 1 # curve-like finite element
            copyto!(csmatout, e1)
        elseif mdim == 2 # surface-like finite element
            n = cross(e1, vec(tangents[:, 2] / norm(tangents[:, 2])))
            e2 = cross(n, e1)
            e2 = e2 / norm(e2)
            csmatout[:, 1] = e1
            csmatout[:, 2] = e2
        end
    end
    return csmatout
end

end
