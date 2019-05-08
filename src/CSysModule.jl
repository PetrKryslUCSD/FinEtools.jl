"""
    CSysModule

Module for management of coordinate systems.
"""
module CSysModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import LinearAlgebra: norm, cross

"""
    CSys{F<:Function}

Type for coordinate system transformations.
"""
struct CSys{F<:Function}
    isconstant::Bool
    isidentity::Bool
    updatebuffer!::F # function to update the coordinate system matrix
    csmat::Array{FFlt, 2} # the coordinate system matrix (buffer); see
end


"""
    CSys(sdim::FInt, mdim::FInt, computecsmat::F) where {F<:Function}

Construct ccoordinate system when the function to compute the
rotation matrix is given.

The function signature:
`update!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)`
where
- `csmatout`= output matrix buffer, 
- `XYZ`= location  in physical coordinates,
- `tangents`= tangent vector matrix, tangents to the parametric coordinate
  curves  in the element,
- `fe_label`= finite element label.
"""
function CSys(sdim::FInt, mdim::FInt, computecsmat::F) where {F<:Function}
    csmat = fill(zero(FFlt), sdim, mdim); # Allocate buffer, in preparation for the first call
    return CSys(false, false, computecsmat, csmat);
end

"""
    CSys(csmat::FFltMat)

Construct coordinate system when the rotation matrix is given.
"""
function CSys(csmat::FFltMat)
    function updatebuffer!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        return csmatout # nothing to be done here, the matrix is already in the buffer
    end
    return CSys(true, false, updatebuffer!, deepcopy(csmat));# fill the buffer with the given matrix
end

"""
    CSys(dim::FInt)

Construct coordinate system when the rotation matrix is the identity.

`dim` = is the space dimension.
"""
function CSys(dim::FInt)
    return CSys([i==j ? one(FFlt) : zero(FFlt) for i=1:dim, j=1:dim]);
end

"""
    CSys(sdim::FInt, mdim::FInt)

Construct coordinate system for isotropic-material used with isoparametric
finite elements.

- `sdim` = number of space dimensions,
- `mdim` = number of manifold dimensions of the finite element in which the
  coordinate system  is being evaluated.

!!! note

If  the coordinate system matrix  should be identity, better use the constructor
for this specific situation, `CSys(dim::FInt)`. That will be much more efficient.

# See also
`gen_iso_csmat`
"""
function CSys(sdim::FInt, mdim::FInt)
    csmat = fill(zero(FFlt), sdim, mdim); # Allocate buffer, prepare for the first call
    function updatebuffer!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
      gen_iso_csmat!(csmatout, XYZ, tangents, fe_label)
      return  csmatout
    end
    return CSys(false, false, updatebuffer!, csmat);
end

"""
    function updatecsmat!(self::CSys, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)

Update the coordinate system orientation matrix.

The  coordinate system matrix is updated based upon the location `XYZ` of the
evaluation point, and possibly on the Jacobian matrix `tangents` within the
element in which the coordinate system matrix is evaluated,  or perhaps on the
label `fe_label` of the finite element.

After this function returns, the coordinate system matrix can be retrieved
from the buffer as `self.csmat`.
"""
function updatecsmat!(self::CSys, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    self.updatebuffer!(self.csmat, XYZ, tangents, fe_label)
    return self.csmat
end

"""
    gen_iso_csmat(XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)

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
function gen_iso_csmat!(csmatout::FFltMat,
    XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    sdim, mdim = size(tangents);
    if sdim == mdim # finite element embedded in space of the same dimension
        copyto!(csmatout, [i==j ? one(FFlt) : zero(FFlt) for i=1:sdim, j=1:sdim]);
    else # lower-dimensional finite element embedded in space of higher dimension
        @assert 0 < mdim < 3
        e1 = tangents[:,1]/norm(tangents[:,1]);
        if mdim == 1 # curve-like finite element
            copyto!(csmatout, e1);
        elseif mdim == 2 # surface-like finite element
            n = cross(e1, vec(tangents[:,2]/norm(tangents[:,2])));
            e2 = cross(n, e1);
            e2 = e2/norm(e2);
            csmatout[:,1] = e1
            csmatout[:,2] = e2
        end
    end
    return csmatout
end

end
