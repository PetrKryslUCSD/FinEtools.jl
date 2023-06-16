"""
    DeforModelRedModule

Module for  model reduction in models of deformation.
"""
module DeforModelRedModule

__precompile__(true)


"""
    AbstractDeforModelRed

Abstract type of model reduction for deformation analysis.
"""
abstract type AbstractDeforModelRed end

"""
    DeforModelRed1D <: AbstractDeforModelRed

Model reduction for 1-D models.
"""
struct DeforModelRed1DStrain <: AbstractDeforModelRed end

"""
    DeforModelRed1D <: AbstractDeforModelRed

Model reduction for 1-D models.
"""
struct DeforModelRed1DStress <: AbstractDeforModelRed end

"""
    DeforModelRed1D

This is alias for `DeforModelRed1DStress`.
"""
DeforModelRed1D = DeforModelRed1DStress

"""
    DeforModelRed2DStrain <: AbstractDeforModelRed

Model reduction for 2-D plane strain models.
"""
struct DeforModelRed2DStrain <: AbstractDeforModelRed end

"""
    DeforModelRed2DStress <: AbstractDeforModelRed

Model reduction for 2-D plane stress models.
"""
struct DeforModelRed2DStress <: AbstractDeforModelRed end

"""
    DeforModelRed2DAxisymm <: AbstractDeforModelRed

Model reduction for 2-D axially symmetric models.
"""
struct DeforModelRed2DAxisymm <: AbstractDeforModelRed end

"""
    DeforModelRed3D <: AbstractDeforModelRed

Model reduction for 3-D models.
"""
struct DeforModelRed3D <: AbstractDeforModelRed end

"""
    nstressstrain(::Type{AbstractDeforModelRed}):: Int

Return the number of stress components involved in the balance equation (and
how many strain components) are there for each model reduction type?
"""
nstressstrain(::Type{DeforModelRed1D}) = 1
nstressstrain(::Type{DeforModelRed2DStrain}) = 3
nstressstrain(::Type{DeforModelRed2DStress}) = 3
nstressstrain(::Type{DeforModelRed2DAxisymm}) = 4
nstressstrain(::Type{DeforModelRed3D}) = 6

"""
    nthermstrain(::Type{DeforModelRed1D})

How many thermal strain components are there for each model reduction type?
"""
nthermstrain(::Type{DeforModelRed1D}) = 1
nthermstrain(::Type{DeforModelRed2DStrain}) = 4
nthermstrain(::Type{DeforModelRed2DStress}) = 3
nthermstrain(::Type{DeforModelRed2DAxisymm}) = 4
nthermstrain(::Type{DeforModelRed3D}) = 6

"""
    stresscomponentmap(::Type{DeforModelRed1D})

Construct a dictionary to map from stress-component symbols to indexes.

# Example
Which component of the stress vector is for instance sigma_x?
Do
```
julia> comp = stresscomponentmap(DeforModelRed1D)
Dict{Symbol,Int64} with 2 entries:
  :xx => 1
  :x  => 1

julia>

julia> comp[:x]
1
```
"""
function stresscomponentmap(::Type{DeforModelRed1D})
    a = Dict{Symbol,Int}()
    a[:x] = 1
    a[:xx] = 1
    return a
end

function stresscomponentmap(::Type{DeforModelRed2DStrain})
    a = Dict{Symbol,Int}()
    a[:x] = 1
    a[:y] = 2
    a[:xy] = 3
    a[:z] = 4
    a[:xx] = 1
    a[:yy] = 2
    a[:zz] = 4
    return a
end

function stresscomponentmap(::Type{DeforModelRed2DStress})
    a = Dict{Symbol,Int}()
    a[:x] = 1
    a[:y] = 2
    a[:xy] = 3
    a[:xx] = 1
    a[:yy] = 2
    return a
end

function stresscomponentmap(::Type{DeforModelRed2DAxisymm})
    a = Dict{Symbol,Int}()
    a[:x] = 1
    a[:y] = 2
    a[:z] = 3
    a[:xy] = 4
    a[:xx] = 1
    a[:yy] = 2
    a[:zz] = 3
    return a
end

function stresscomponentmap(::Type{DeforModelRed3D})
    a = Dict{Symbol,Int}()
    a[:x] = 1
    a[:y] = 2
    a[:z] = 3
    a[:xy] = 4
    a[:xz] = 5
    a[:yz] = 6
    a[:xx] = 1
    a[:yy] = 2
    a[:zz] = 3
    return a
end

"""
    blmat!(
        MR::Type{DeforModelRed1D},
        B::AbstractMatrix{T},
        N::AbstractMatrix{T},
        gradN::AbstractMatrix{T},
        c::AbstractMatrix{T},
        Rm::AbstractMatrix{T},
    ) where {T<:Number}

Compute the strain-displacement matrix for a one-manifold element.

Compute the linear, displacement independent, strain-displacement matrix for a
one-manifold element.   *The input displacements are in the global Cartesian
coordinate system, the output strains are in the material coordinate system.*

# Arguments
- `N` =matrix of basis function values
- `gradN` =matrix of basis function gradients with respect to the Cartesian
  coordinates in the directions of the material orientation
- `c` =array of spatial coordinates of the evaluation point in the global
  Cartesian coordinates.
- `Rm` =orthogonal matrix with the unit basis vectors of the local material
  orientation coordinate system as columns. `size(Rm)= (ndim,mdim)`, where
  `ndim` = number of spatial dimensions of the embedding space (here `ndim <=
  3`), and `mdim` = number of manifold dimensions (here `mdim == 1`).

# Output
- `B` = strain-displacement matrix, where  `size(B) =
  (nstressstrain,nnodes*dim)`; here `nstressstrain`= number of strains, `dim`
  = number of spatial dimensions of the embedding space, and `nnodes` = number
  of finite element nodes on the element. The strain components are ordered as
  shown  in `stresscomponentmap`. The matrix is passed in as a buffer, set to
  zero,  and filled in  with the nonzero components.  It is also returned for
  convenience.
"""
function blmat!(
    MR::Type{DeforModelRed1D},
    B::AbstractMatrix{T},
    N::AbstractMatrix{T},
    gradN::AbstractMatrix{T},
    c::AbstractMatrix{T},
    Rm::AbstractMatrix{T},
) where {T<:Number}
    nnodes = size(gradN, 1)
    @assert (dim = size(c, 2)) >= 1
    @assert size(B) == (nstressstrain(MR), dim * nnodes)
    fill!(B, zero(T))
    @inbounds for i in 1:nnodes
        for j in 1:dim
            k = dim * (i - 1)
            B[1, k+j] = gradN[i, 1] * Rm[j, 1]
        end
    end
    return B
end

"""
    blmat!(
        MR::Type{DeforModelRed2DStrain},
        B::AbstractMatrix{T},
        N::AbstractMatrix{T},
        gradN::AbstractMatrix{T},
        c::AbstractMatrix{T},
        Rm::AbstractMatrix{T},
    ) where {T<:Number}

Compute the strain-displacement matrix for a two-manifold element for
plane-strain models.

Compute the linear, displacement independent, strain-displacement matrix for a
two-manifold element in plane-strain setting.   *The input displacements are
in the global Cartesian coordinate system, the output strains are in the
material coordinate system.*

# Arguments
- `N` =matrix of basis function values
- `gradN` =matrix of basis function gradients with respect to the Cartesian
  coordinates in the directions of the material orientation
- `c` =array of spatial coordinates of the evaluation point in the global
  Cartesian coordinates.
- `Rm` =orthogonal matrix with the unit basis vectors of the local material
  orientation coordinate system as columns. `size(Rm)= (ndim,mdim)`, where
  `ndim` = number of spatial dimensions of the embedding space (here `ndim <=
  3`), and `mdim` = number of manifold dimensions (here `mdim == 2`).

# Output
- `B` = strain-displacement matrix, where  `size(B) =
  (nstressstrain,nnodes*dim)`; here `nstressstrain`= number of strains, `dim`
  = number of spatial dimensions of the embedding space, and `nnodes` = number
  of finite element nodes on the element. The strain components are ordered as
  shown  in `stresscomponentmap`. The matrix is passed in as a buffer, set to
  zero,  and filled in  with the nonzero components.  It is also returned for
  convenience.
"""
function blmat!(
    MR::Type{DeforModelRed2DStrain},
    B::AbstractMatrix{T},
    N::AbstractMatrix{T},
    gradN::AbstractMatrix{T},
    c::AbstractMatrix{T},
    Rm::AbstractMatrix{T},
) where {T<:Number}
    nnodes = size(gradN, 1)
    @assert (dim = size(c, 2)) == 2
    @assert size(B) == (nstressstrain(MR), dim * nnodes)
    fill!(B, zero(T))
    @inbounds for i in 1:nnodes
        k = dim * (i - 1)
        @inbounds for j in 1:dim
            B[1, k+j] = gradN[i, 1] * Rm[j, 1]
            B[2, k+j] = gradN[i, 2] * Rm[j, 2]
            B[3, k+j] = gradN[i, 2] * Rm[j, 1] + gradN[i, 1] * Rm[j, 2]
        end
    end
    return B
end

"""
    blmat!(
        MR::Type{DeforModelRed2DStress},
        B::AbstractMatrix{T},
        N::AbstractMatrix{T},
        gradN::AbstractMatrix{T},
        c::AbstractMatrix{T},
        Rm::AbstractMatrix{T},
    ) where {T<:Number}

Compute the strain-displacement matrix for a two-manifold element for
plane-stress models.

Compute the linear, displacement independent, strain-displacement matrix for a
two-manifold element in plane-stress setting.   *The input displacements
are in the global Cartesian coordinate system, the output strains are in the
material coordinate system.*

# Arguments
- `N` =matrix of basis function values
- `gradN` =matrix of basis function gradients with respect to the Cartesian
  coordinates in the directions of the material orientation
- `c` =array of spatial coordinates of the evaluation point in the global
  Cartesian coordinates.
- `Rm` =orthogonal matrix with the unit basis vectors of the local material
  orientation coordinate system as columns. `size(Rm)= (ndim,mdim)`, where
  `ndim` = number of spatial dimensions of the embedding space (here `ndim <=
  3`), and `mdim` = number of manifold dimensions (here `mdim == 2`).

# Output
- `B` = strain-displacement matrix, where  `size(B) =
  (nstressstrain,nnodes*dim)`; here `nstressstrain`= number of strains, `dim`
  = number of spatial dimensions of the embedding space, and `nnodes` = number
  of finite element nodes on the element. The strain components are ordered as
  shown  in `stresscomponentmap`. The matrix is passed in as a buffer, set to
  zero,  and filled in  with the nonzero components.  It is also returned for
  convenience.
"""
function blmat!(
    MR::Type{DeforModelRed2DStress},
    B::AbstractMatrix{T},
    N::AbstractMatrix{T},
    gradN::AbstractMatrix{T},
    c::AbstractMatrix{T},
    Rm::AbstractMatrix{T},
) where {T<:Number}
    nnodes = size(gradN, 1)
    @assert (dim = size(c, 2)) >= 2
    @assert size(B) == (nstressstrain(MR), dim * nnodes)
    fill!(B, zero(T))
    @inbounds for i in 1:nnodes
        k = dim * (i - 1)
        for j in 1:dim
            B[1, k+j] = gradN[i, 1] * Rm[j, 1]
            B[2, k+j] = gradN[i, 2] * Rm[j, 2]
            B[3, k+j] = gradN[i, 2] * Rm[j, 1] + gradN[i, 1] * Rm[j, 2]
        end
    end
    return B
end

"""
    blmat!(
        MR::Type{DeforModelRed2DAxisymm},
        B::AbstractMatrix{T},
        N::AbstractMatrix{T},
        gradN::AbstractMatrix{T},
        c::AbstractMatrix{T},
        Rm::AbstractMatrix{T},
    ) where {T<:Number}

Compute the strain-displacement matrix for a two-manifold element for
axially symmetric models.

Compute the linear, displacement independent, strain-displacement matrix for a
two-manifold element in axially-symmetric setting.   *The input displacements
are in the global Cartesian coordinate system, the output strains are in the
material coordinate system.*

# Arguments
- `N` =matrix of basis function values
- `gradN` =matrix of basis function gradients with respect to the Cartesian
  coordinates in the directions of the material orientation
- `c` =array of spatial coordinates of the evaluation point in the global
  Cartesian coordinates.
- `Rm` =orthogonal matrix with the unit basis vectors of the local material
  orientation coordinate system as columns. `size(Rm)= (ndim,mdim)`, where
  `ndim` = number of spatial dimensions of the embedding space (here `ndim <=
  3`), and `mdim` = number of manifold dimensions (here `mdim == 2`).

# Output
- `B` = strain-displacement matrix, where  `size(B) =
  (nstressstrain,nnodes*dim)`; here `nstressstrain`= number of strains, `dim`
  = number of spatial dimensions of the embedding space, and `nnodes` = number
  of finite element nodes on the element. The strain components are ordered as
  shown  in `stresscomponentmap`. The matrix is passed in as a buffer, set to
  zero,  and filled in  with the nonzero components.  It is also returned for
  convenience.
"""
function blmat!(
    MR::Type{DeforModelRed2DAxisymm},
    B::AbstractMatrix{T},
    N::AbstractMatrix{T},
    gradN::AbstractMatrix{T},
    c::AbstractMatrix{T},
    Rm::AbstractMatrix{T},
) where {T<:Number}
    nnodes = size(gradN, 1)
    @assert (dim = size(c, 2)) == 2
    @assert size(B) == (nstressstrain(MR), dim * nnodes)
    r = c[1] # distance from the axis of symmetry
    # if practically zero, clamped it to machine epsilon because we need to divide with it
    if r <= eps(1.0)
        r = eps(1.0)
    end
    fill!(B, zero(T))
    @inbounds for i in 1:nnodes
        k = dim * (i - 1)
        for j in 1:dim
            B[1, k+j] = gradN[i, 1] * Rm[j, 1]
            B[2, k+j] = gradN[i, 2] * Rm[j, 2]
            B[3, k+j] = N[i] / r * Rm[j, 1]
            B[4, k+j] = gradN[i, 2] * Rm[j, 1] + gradN[i, 1] * Rm[j, 2]
        end
    end
    return B
end

"""
    blmat!(
        MR::Type{DeforModelRed3D},
        B::AbstractMatrix{T},
        N::AbstractMatrix{T},
        gradN::AbstractMatrix{T},
        c::AbstractMatrix{T},
        Rm::AbstractMatrix{T},
    ) where {T<:Number}

Compute the strain-displacement matrix for a three-manifold element.

Compute the linear, displacement independent, strain-displacement matrix
for a three-manifold element.   *The input displacements are in the global
Cartesian coordinate system, the output strains are in the material coordinate
system.*

# Arguments
- `N` =matrix of basis function values
- `gradN` =matrix of basis function gradients with respect to the Cartesian
  coordinates in the directions of the material orientation
- `c` =array of spatial coordinates of the evaluation point in the global
  Cartesian coordinates.
- `Rm` =orthogonal matrix with the unit basis vectors of the local material
  orientation coordinate system as columns. `size(Rm)= (ndim,mdim)`, where
  `ndim` = number of spatial dimensions of the embedding space (here `ndim ==
  3`), and `mdim` = number of manifold dimensions (here `mdim == 3`).

# Output
- `B` = strain-displacement matrix, where  `size(B) =
  (nstressstrain,nnodes*dim)`; here `nstressstrain`= number of strains, `dim`
  = number of spatial dimensions of the embedding space, and `nnodes` = number
  of finite element nodes on the element. The strain components are ordered as
  shown  in `stresscomponentmap`. The matrix is passed in as a buffer, set to
  zero,  and filled in  with the nonzero components.  It is also returned for
  convenience.
"""
function blmat!(
    MR::Type{DeforModelRed3D},
    B::AbstractMatrix{T},
    N::AbstractMatrix{T},
    gradN::AbstractMatrix{T},
    c::AbstractMatrix{T},
    Rm::AbstractMatrix{T},
) where {T<:Number}
    nnodes = size(gradN, 1)
    @assert (dim = size(c, 2)) == 3
    @assert size(B) == (nstressstrain(MR), dim * nnodes)
    fill!(B, zero(T))
    @inbounds for i in 1:nnodes
        k = dim * (i - 1)
        for j in 1:dim
            kj = k + j
            B[1, kj] = gradN[i, 1] * Rm[j, 1]
            B[2, kj] = gradN[i, 2] * Rm[j, 2]
            B[3, kj] = gradN[i, 3] * Rm[j, 3]
            B[4, kj] = gradN[i, 2] * Rm[j, 1] + gradN[i, 1] * Rm[j, 2]
            B[5, kj] = gradN[i, 3] * Rm[j, 1] + gradN[i, 1] * Rm[j, 3]
            B[6, kj] = gradN[i, 3] * Rm[j, 2] + gradN[i, 2] * Rm[j, 3]
        end
    end
    return B
end

end # module


# """
#     divmat(MR::Type{DeforModelRed3D}, N::AbstractMatrix{T}, gradN::AbstractMatrix{T}, c::AbstractMatrix{T})

# Compute the displacement divergence matrix for a three-manifold element.


# # Arguments
# - `N` =matrix of basis function values
# - `gradN` =matrix of basis function gradients with respect to the Cartesian
#   coordinates in the directions of the material orientation
# - `c` =array of spatial coordinates of the evaluation point in the global
#   Cartesian coordinates.
# - `Rm` =orthogonal matrix with the unit basis vectors of the local material
#   orientation coordinate system as columns. `size(Rm)= (ndim,mdim)`, where
#   `ndim` = number of spatial dimensions of the embedding space (here `ndim ==
#   3`), and `mdim` = number of manifold dimensions (here `mdim == 3`).

# # Output
# - `divmat` = displacement divergence matrix, where  `size(divmat) =
#   (1,nnodes*dim)`; here `dim` = Number of spatial dimensions of the embedding
#   space, and `nnodes` = number of finite element nodes on the element.  The matrix
#   is passed in as a buffer, set to zero,  and filled in  with the nonzero
#   components.  It is also returned for convenience.
# """
# function divmat(MR::Type{DeforModelRed3D}, N::AbstractMatrix{T}, gradN::AbstractMatrix{T}, c::AbstractMatrix{T})
#     nnodes = size(gradN, 1)
#     @assert (dim = size(c, 2)) == 3
#     divm = fill(0.0, 1, dim * nnodes)
#     for j = 1:dim
#         for i = 1:nnodes
#             k = dim * (i - 1)
#             divm[1, k+j] = gradN[i, j]
#         end
#     end
#     return divm::AbstractMatrix{T}
# end

# """
#     vgradmat(MR::Type{DeforModelRed3D}, N::AbstractMatrix{T}, gradN::AbstractMatrix{T}, c::AbstractMatrix{T})

# Compute the matrix of displacement gradient in vector form for a three-manifold element.


# # Arguments
# - `N` =matrix of basis function values
# - `gradN` =matrix of basis function gradients with respect to the Cartesian
#   coordinates in the directions of the material orientation
# - `c` =array of spatial coordinates of the evaluation point in the global
#   Cartesian coordinates.
# - `Rm` =orthogonal matrix with the unit basis vectors of the local material
#   orientation coordinate system as columns. `size(Rm)= (ndim,mdim)`, where
#   `ndim` = number of spatial dimensions of the embedding space (here `ndim ==
#   3`), and `mdim` = number of manifold dimensions (here `mdim == 3`).

# # Output
# - `vgradm` = strain-displacement matrix, where  `size(divmat) =
#   (1,nnodes*dim)`; here `dim` = Number of spatial dimensions of the embedding
#   space, and `nnodes` = number of finite element nodes on the element.  The matrix
#   is passed in as a buffer, set to zero,  and filled in  with the nonzero
#   components.  It is also returned for convenience.
# """
# function vgradmat(MR::Type{DeforModelRed3D}, N::AbstractMatrix{T}, gradN::AbstractMatrix{T}, c::AbstractMatrix{T})
#     nnodes = size(gradN, 1)
#     @assert (dim = size(c, 2)) == 3
#     vgradm = fill(0.0, dim * dim, dim * nnodes)
#     for i = 1:dim
#         vgradm[dim*(i-1)+1:dim*i, i:dim:nnodes*dim-dim+i] .= transpose(gradN)
#     end
#     return vgradm::AbstractMatrix{T}
# end

