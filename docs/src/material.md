[Table of contents](https://petrkryslucsd.github.io/FinEtools.jl)

# Material and Material Orientation

The material response  is described in  material-point-attached coordinate system. These coordinate systems  are Cartesian, and the material coordinate system is typically chosen to make  the response particularly simple.  So for orthotropic or transversely isotropic materials the axes would be aligned with the axes of orthotropy.

The type `CSys` (module `CSysModule`) is the updater of the material coordinate system matrix. The object is equipped with a callback to store the current orientation matrix. For instance: the coordinate system for an orthotropic material wound around a cylinder could be described in the coordinate system `CSys(3, 3, updatecs!)`, where the callback is defined as

```julia
function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    csmatout[:, 2] = [0.0 0.0 1.0]
    csmatout[:, 3] = XYZ
    csmatout[3, 3] = 0.0
    csmatout[:, 3] = csmatout[:, 3]/norm(csmatout[:, 3])
    csmatout[:, 1] = cross(csmatout[:, 2], csmatout[:, 3])
end
```

## Materials for linear deformation analysis

The module `MatDeforModule` provides functions to convert between vector and matrix (tensor) representations of stress and strain. Further, functions to rotate stress and strain between different coordinate systems (based upon the model-reduction type `DeforModelRed`) are provided.

Currently  there are material types for isotropic and orthotropic linear elastic materials. The user may add  additional material types by deriving from `MatDefor` and equipping them with three methods: (1) compute the tangent moduli, (2) update the material state, (3) compute the thermal strain.

For full generality, material types  should implement these methods for fully three-dimensional, plane strain and plane stress, 2D axially symmetric, and one-dimensional deformation models.
