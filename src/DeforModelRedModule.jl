module DeforModelRedModule

export DeforModelRed, DeforModelRed1D, DeforModelRed2DStrain,
       DeforModelRed2DStress, DeforModelRed2DAxisymm, DeforModelRed3D
export nstsstn, nthstn, noutsts
export stresscomponentmap
export Blmat!

using FinEtools.FTypesModule

abstract type DeforModelRed end;
struct DeforModelRed1D <: DeforModelRed
end
struct DeforModelRed2DStrain <: DeforModelRed
end
struct DeforModelRed2DStress <: DeforModelRed
end
struct DeforModelRed2DAxisymm <: DeforModelRed
end
struct DeforModelRed3D <: DeforModelRed
end

"""
    nstsstn(::Type{DeforModelRed1D})::FInt

How many stress components involved in the balance equation (and how many
strain components) are there for each model reduction type?
"""
nstsstn(::Type{DeforModelRed1D})::FInt = 1;
nstsstn(::Type{DeforModelRed2DStrain})::FInt = 3;
nstsstn(::Type{DeforModelRed2DStress})::FInt = 3;
nstsstn(::Type{DeforModelRed2DAxisymm})::FInt = 4;
nstsstn(::Type{DeforModelRed3D})::FInt = 6;

"""
    nthstn(::Type{DeforModelRed1D})::FInt

How many thermal strain components are there for each model reduction type?
"""
nthstn(::Type{DeforModelRed1D})::FInt = 1;
nthstn(::Type{DeforModelRed2DStrain})::FInt = 4;
nthstn(::Type{DeforModelRed2DStress})::FInt = 3;
nthstn(::Type{DeforModelRed2DAxisymm})::FInt = 4;
nthstn(::Type{DeforModelRed3D})::FInt = 6;
 
"""
    stresscomponentmap(::Type{DeforModelRed1D})

Construct a dictionary to map from stress-component symbols to indexes.
"""
function stresscomponentmap(::Type{DeforModelRed1D})
    a=Dict{Symbol,FInt}();
    a[:x]=1;
    a[:xx]=1;
    return a
end

function stresscomponentmap(::Type{DeforModelRed2DStrain})
    a=Dict{Symbol,FInt}();
    a[:x]=1;  a[:y]=2;  a[:xy]=3;  a[:z]=4;
    a[:xx]=1;  a[:yy]=2; a[:zz]=4;
    return a
end

function stresscomponentmap(::Type{DeforModelRed2DStress})
    a=Dict{Symbol,FInt}();
    a[:x]=1;  a[:y]=2;  a[:xy]=3;
    a[:xx]=1;  a[:yy]=2;
    return a
end

function stresscomponentmap(::Type{DeforModelRed2DAxisymm})
    a=Dict{Symbol,FInt}();
    a[:x]=1;  a[:y]=2;  a[:z]=3;  a[:xy]=4;
    a[:xx]=1;  a[:yy]=2;  a[:zz]=3;
    return a
end

function stresscomponentmap(::Type{DeforModelRed3D})
    a=Dict{Symbol,FInt}();
    a[:x]=1;  a[:y]=2;  a[:z]=3;  a[:xy]=4; a[:xz]=5; a[:yz]=6;
    a[:xx]=1;  a[:yy]=2;  a[:zz]=3;
    return a
end

"""
    Blmat!(MR::Type{DeforModelRed1D}, B::FFltMat,
      N::FFltMat, gradN::FFltMat, c::FFltMat, Rm::FFltMat)

Compute the strain-displacement matrix for a one-manifold element.

Compute the linear, displacement independent, strain-displacement matrix for a
one-manifold element.   *The input displacements are in the global Cartesian
coordinate system, the output strains are in the material coordinate system.*

`N` =matrix of basis function values
`gradN` =matrix of basis function gradients with respect to the
  Cartesian coordinates in the directions of the material orientation
`c` =array of spatial coordinates of the evaluation point
  in the global Cartesian coordinates.
`Rm` =orthogonal matrix with the unit basis vectors of the local
  material orientation coordinate system as columns.
  `size(Rm)= [ndim,mdim]`, where `ndim` = number of spatial dimensions
  of the embedding space (here `ndim >= 3`), aand `mdim` = number of manifold
  dimensions (here `mdim == 1`).

Output:
`B` = strain-displacement matrix, where  `size(B) = [nstsstn,nnodes*dim]`;
  here nstsstn= number of strains, dim = Number of spatial
  dimensions of the embedding space, and nnodes = number of
  finite element nodes on the element. The strain components are ordered as
  shown  in `stresscomponentmap`. The matrix is passed in as a buffer,
  set to zero,  and filled in  with the nonzero components.  It is also returned
for convenience.
"""
function Blmat!(MR::Type{DeforModelRed1D}, B::FFltMat,
  N::FFltMat, gradN::FFltMat, c::FFltMat, Rm::FFltMat)
  nnodes = size(gradN,1);
  @assert (dim = size(c,2)) >= 1;
  @assert size(B) == (nstsstn(MR), dim*nnodes)
  fill!(B,0.0);
  @inbounds for i= 1:nnodes
    @inbounds for j=1:dim
      k=dim*(i-1);
      B[1,k+j]=  gradN[i,1] *Rm[j,1];
    end
  end
  return B::FFltMat;
end

"""
    Blmat!(MR::Type{DeforModelRed2DStrain}, B::FFltMat,
      N::FFltMat, gradN::FFltMat, c::FFltMat, Rm::FFltMat)

Compute the strain-displacement matrix for a two-manifold element for
plane-strain models.

Compute the linear, displacement independent, strain-displacement matrix for a
two-manifold element in plane-strain setting.   *The input displacements
are in the global Cartesian coordinate system, the output strains are in the
material coordinate system.*

`N` =matrix of basis function values
`gradN` =matrix of basis function gradients with respect to the
  Cartesian coordinates in the directions of the material orientation
`c` =array of spatial coordinates of the evaluation point
  in the global Cartesian coordinates.
`Rm` =orthogonal matrix with the unit basis vectors of the local
  material orientation coordinate system as columns.
  `size(Rm)= [ndim,mdim]`, where `ndim` = number of spatial dimensions
  of the embedding space (here `ndim == 2`), aand `mdim` = number of manifold
  dimensions (here `mdim == 2`).

Output:
`B` = strain-displacement matrix, where  `size(B) = [nstsstn,nnodes*dim]`;
  here nstsstn= number of strains, dim = Number of spatial
  dimensions of the embedding space, and nnodes = number of
  finite element nodes on the element. The strain components are ordered as
  shown  in `stresscomponentmap`. The matrix is passed in as a buffer,
  set to zero,  and filled in  with the nonzero components.  It is also returned
for convenience.
"""
function Blmat!(MR::Type{DeforModelRed2DStrain}, B::FFltMat,
  N::FFltMat, gradN::FFltMat, c::FFltMat, Rm::FFltMat)
  nnodes = size(gradN,1);
  @assert (dim = size(c,2)) == 2;
  @assert size(B) == (nstsstn(MR), dim*nnodes)
  fill!(B,0.0);
  @inbounds for i= 1:nnodes
    k=dim*(i-1);
    @inbounds for j=1:dim
      B[1,k+j]=gradN[i,1]*Rm[j,1]
      B[2,k+j]=gradN[i,2]*Rm[j,2]
      B[3,k+j]=gradN[i,2]*Rm[j,1]+gradN[i,1]*Rm[j,2]
    end
  end
  return B::FFltMat;
end

"""
    Blmat!(MR::Type{DeforModelRed2DStress}, B::FFltMat,
      N::FFltMat, gradN::FFltMat, c::FFltMat, Rm::FFltMat)

Compute the strain-displacement matrix for a two-manifold element for
plane-stress models.

Compute the linear, displacement independent, strain-displacement matrix for a
two-manifold element in plane-stress setting.   *The input displacements
are in the global Cartesian coordinate system, the output strains are in the
material coordinate system.*

`N` =matrix of basis function values
`gradN` =matrix of basis function gradients with respect to the
  Cartesian coordinates in the directions of the material orientation
`c` =array of spatial coordinates of the evaluation point
  in the global Cartesian coordinates.
`Rm` =orthogonal matrix with the unit basis vectors of the local
  material orientation coordinate system as columns.
  `size(Rm)= [ndim,mdim]`, where `ndim` = number of spatial dimensions
  of the embedding space (here `ndim >= 3`), aand `mdim` = number of manifold
  dimensions (here `mdim=2`).

Output:
`B` = strain-displacement matrix, where  `size(B) = [nstsstn,nnodes*dim]`;
  here nstsstn= number of strains, dim = Number of spatial
  dimensions of the embedding space, and nnodes = number of
  finite element nodes on the element. The strain components are ordered as
  shown  in `stresscomponentmap`. The matrix is passed in as a buffer,
  set to zero,  and filled in  with the nonzero components.  It is also returned
for convenience.
"""
function Blmat!(MR::Type{DeforModelRed2DStress}, B::FFltMat,
  N::FFltMat, gradN::FFltMat, c::FFltMat, Rm::FFltMat)
  nnodes = size(gradN,1);
  @assert (dim = size(c,2)) >= 2;
  @assert size(B) == (nstsstn(MR), dim*nnodes)
  fill!(B,0.0);
  @inbounds for i= 1:nnodes
    k=dim*(i-1);
    @inbounds for j=1:dim
      B[1,k+j]=gradN[i,1]*Rm[j,1]
      B[2,k+j]=gradN[i,2]*Rm[j,2]
      B[3,k+j]=gradN[i,2]*Rm[j,1]+gradN[i,1]*Rm[j,2]
    end
  end
  return B::FFltMat;
end

"""
    Blmat!(MR::Type{DeforModelRed2DAxisymm}, B::FFltMat,
      N::FFltMat, gradN::FFltMat, c::FFltMat, Rm::FFltMat)

Compute the strain-displacement matrix for a two-manifold element for
axially symmetric models.

Compute the linear, displacement independent, strain-displacement matrix for a
two-manifold element in axially-symmetric setting.   *The input displacements
are in the global Cartesian coordinate system, the output strains are in the
material coordinate system.*

`N` =matrix of basis function values
`gradN` =matrix of basis function gradients with respect to the
  Cartesian coordinates in the directions of the material orientation
`c` =array of spatial coordinates of the evaluation point
  in the global Cartesian coordinates; `c[1]` is the distance from the axis of
  symmetry.
`Rm` =orthogonal matrix with the unit basis vectors of the local
  material orientation coordinate system as columns.
  `size(Rm)= [ndim,ndim]`, where `ndim` = number of spatial dimensions
  of the embedding space (here `ndim=3`).

Output:
`B` = strain-displacement matrix, where  `size(B) = [nstsstn,nnodes*dim]`;
  here nstsstn= number of strains, dim = Number of spatial
  dimensions of the embedding space, and nnodes = number of
  finite element nodes on the element. The strain components are ordered as
  shown  in `stresscomponentmap`. The matrix is passed in as a buffer,
  set to zero,  and filled in  with the nonzero components.  It is also returned
for convenience.
"""
function Blmat!(MR::Type{DeforModelRed2DAxisymm}, B::FFltMat,
  N::FFltMat, gradN::FFltMat, c::FFltMat, Rm::FFltMat)
  nnodes = size(gradN,1);
  @assert (dim = size(c,2)) == 2;
  @assert size(B) == (nstsstn(MR), dim*nnodes)
  r = c[1]; # distance from the axis of symmetry
  # if practically zero, clamped it to machine epsilon because we need to divide with it
  if r <= eps(1.0)
    r = eps(1.0);
  end
  fill!(B,0.0);
  @inbounds for i= 1:nnodes
    k=dim*(i-1);
    @inbounds for j=1:dim
      B[1,k+j]=gradN[i,1]*Rm[j,1]
      B[2,k+j]=gradN[i,2]*Rm[j,2]
      B[3,k+j]=N[i]/r*Rm[j,1]
      B[4,k+j]=gradN[i,2]*Rm[j,1]+gradN[i,1]*Rm[j,2]
    end
  end
  return B::FFltMat;
end

"""
    Blmat!(MR::Type{DeforModelRed3D}, B::FFltMat,
      N::FFltMat, gradN::FFltMat, c::FFltMat, Rm::FFltMat)

Compute the strain-displacement matrix for a three-manifold element.

Compute the linear, displacement independent, strain-displacement matrix
for a three-manifold element.   *The input displacements are in the global
Cartesian coordinate system, the output strains are in the material coordinate
system.*

`N` =matrix of basis function values
`gradN` =matrix of basis function gradients with respect to the
  Cartesian coordinates in the directions of the material orientation
`c` =array of spatial coordinates of the evaluation point
  in the global Cartesian coordinates
`Rm` =orthogonal matrix with the unit basis vectors of the local
  material orientation coordinate system as columns.
  `size(Rm)= [ndim,ndim]`, where `ndim` = number of spatial dimensions
  of the embedding space (here `ndim=3`).

Output:
`B` = strain-displacement matrix, where  `size(B) = [nstsstn,nnodes*dim]`;
  here nstsstn= number of strains, dim = Number of spatial
  dimensions of the embedding space, and nnodes = number of
  finite element nodes on the element. The strain components are ordered as
  shown  in `stresscomponentmap`. The matrix is passed in as a buffer,
  set to zero,  and filled in  with the nonzero components.  It is also returned
for convenience.
"""
function Blmat!(MR::Type{DeforModelRed3D}, B::FFltMat,
  N::FFltMat, gradN::FFltMat, c::FFltMat, Rm::FFltMat)
  nnodes = size(gradN,1);
  @assert (dim = size(c,2)) == 3;
  @assert size(B) == (nstsstn(MR), dim*nnodes)
  fill!(B,0.0);
  @inbounds for i= 1:nnodes
    k=dim*(i-1);
    @inbounds for j=1:dim
      kj=k+j
      B[1,kj]=gradN[i,1]*Rm[j,1]
      B[2,kj]=gradN[i,2]*Rm[j,2]
      B[3,kj]=gradN[i,3]*Rm[j,3]
      B[4,kj]=gradN[i,2]*Rm[j,1]+gradN[i,1]*Rm[j,2]
      B[5,kj]=gradN[i,3]*Rm[j,1]+gradN[i,1]*Rm[j,3]
      B[6,kj]=gradN[i,3]*Rm[j,2]+gradN[i,2]*Rm[j,3]
    end
  end
  return B::FFltMat;
end


end
