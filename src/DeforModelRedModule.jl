module DeforModelRedModule

using FinEtools.FTypesModule


abstract type DeforModelRed end;
immutable DeforModelRed1D <: DeforModelRed
end
immutable DeforModelRed2DStrain <: DeforModelRed
end
immutable DeforModelRed2DStress <: DeforModelRed
end
immutable DeforModelRed2DAxisymm <: DeforModelRed
end
immutable DeforModelRed3D <: DeforModelRed
end
export DeforModelRed
export DeforModelRed1D
export DeforModelRed2DStrain
export DeforModelRed2DStress
export DeforModelRed2DAxisymm
export DeforModelRed3D

# How many stress components involved in the balance equation are there for each
# model reduction type?
nstresses(::Type{DeforModelRed1D})::FInt = 1;
export nstresses
nstresses(::Type{DeforModelRed2DStrain})::FInt = 3;
export nstresses
nstresses(::Type{DeforModelRed2DStress})::FInt = 3;
export nstresses
nstresses(::Type{DeforModelRed2DAxisymm})::FInt = 4;
export nstresses
nstresses(::Type{DeforModelRed3D})::FInt = 6;
export nstresses
# How many thermal strain components are there for each model reduction type?
nthstrains(::Type{DeforModelRed1D})::FInt = 1;
export nthstrains
nthstrains(::Type{DeforModelRed2DStrain})::FInt = 4;
export nthstrains
nthstrains(::Type{DeforModelRed2DStress})::FInt = 3;
export nthstrains
nthstrains(::Type{DeforModelRed2DAxisymm})::FInt = 4;
export nthstrains
nthstrains(::Type{DeforModelRed3D})::FInt = 6;
export nthstrains

function stresscomponentmap(::Type{DeforModelRed1D})
    a=Dict{Symbol,Int64}(); a[:x]=1;
    a[:xx]=1;
    return a
end
export stresscomponentmap

function stresscomponentmap(::Type{DeforModelRed2DStrain})
    a=Dict{Symbol,Int64}(); a[:x]=1;  a[:y]=2;  a[:xy]=3;  a[:z]=4;
    a[:xx]=1;  a[:yy]=2; a[:zz]=4;
    return a
end
export stresscomponentmap

function stresscomponentmap(::Type{DeforModelRed2DStress})
    a=Dict{Symbol,Int64}(); a[:x]=1;  a[:y]=2;  a[:xy]=3;
    a[:xx]=1;  a[:yy]=2;
    return a
end
export stresscomponentmap

function stresscomponentmap(::Type{DeforModelRed2DAxisymm})
    a=Dict{Symbol,Int64}(); a[:x]=1;  a[:y]=2;  a[:z]=3;  a[:xy]=4;
    a[:xx]=1;  a[:yy]=2;  a[:zz]=3;
    return a
end
export stresscomponentmap

function stresscomponentmap(::Type{DeforModelRed3D})
    a=Dict{Symbol,Int64}(); a[:x]=1;  a[:y]=2;  a[:z]=3;  a[:xy]=4; a[:xz]=5; a[:yz]=6;
    a[:xx]=1;  a[:yy]=2;  a[:zz]=3;
    return a
end
export stresscomponentmap

function Blmat!(::Type{DeforModelRed1D},B::FFltMat,N::FFltMat,gradN::FFltMat,c::FFltMat,Rm::FFltMat)
    # Compute the strain-displacement matrix for a one-manifold element.
    #
    # B = Blmat1(self,N,gradN,c,Rm)
    #
    # Compute the linear, displacement independent, strain-displacement matrix
    # for a one-manifold element.
    #
    # Arguments
    #      self=finite element model
    #      N= matrix of basis function values
    #      gradN=matrix of basis function gradients with respect to the
    #         Cartesian coordinates in the directions of the material orientation
    #      c=array of spatial coordinates of the evaluation point
    #         in the global Cartesian coordinates
    #      Rm=orthogonal matrix with the unit basis vectors of the local
    #         material orientation coordinate system as columns;
    #         supply Rm as empty when there is no need for global-to-local
    #         transformation (for instance for isotropic materials).
    #         size(Rm)= [ndim,2], where ndim = number of spatial dimensions
    #         of the embedding space, greater than or equal to  2.
    # Output:
    #      B = strain-displacement matrix, size (B) = [1,nnodes*ndim],
    #         where  nnodes = number of
    #         finite element nodes on the element.
    #
    nfn= size(gradN,1);
    dim =size(c,2);
    @assert size(B)==(1,dim*nfn)
    fill!(B,0.0);
    @inbounds for i= 1:nfn
        @inbounds for j=1:dim
            k=dim*(i-1);
            B[1,k+j]=  gradN[i,1] *Rm[j,1];
        end
    end
    return B::FFltMat;
end
export Blmat!

function Blmat!(::Type{DeforModelRed2DStrain},
                B::FFltMat,N::FFltMat,gradN::FFltMat,c::FFltMat,Rm::FFltMat)
    # Compute the strain-displacement matrix for a two-manifold element.
    #
    # B = Blmat2(self,N,gradN,c,Rm)
    #
    # Compute the linear, displacement independent, strain-displacement matrix
    # for a two-manifold element.
    #
    # Arguments
    #      self=finite element model
    #      N= matrix of basis function values
    #      gradN=matrix of basis function gradients with respect to the
    #         Cartesian coordinates in the directions of the material orientation
    #      c=array of spatial coordinates of the evaluation point
    #         in the global Cartesian coordinates
    #      Rm=orthogonal matrix with the unit basis vectors of the local
    #         material orientation coordinate system as columns;
    #         supply Rm as empty when there is no need for global-to-local
    #         transformation (for instance for isotropic materials).
    #         size(Rm)= [ndim,2], where ndim = number of spatial dimensions
    #         of the embedding space, greater than or equal to  2.
    # Output:
    #      B = strain-displacement matrix, size (B) = [3,nnodes*2],
    #         where  nnodes = number of
    #         finite element nodes on the element.
    #
    nfn= size(gradN,1);
    dim =size(c,2);
    @assert size(B)==(3,dim*nfn)
    fill!(B,0.0);
    @inbounds for i= 1:nfn
        k=dim*(i-1);
        @inbounds for j=1:dim
            B[1,k+j]=gradN[i,1]*Rm[j,1]
            B[2,k+j]=gradN[i,2]*Rm[j,2]
            B[3,k+j]=gradN[i,2]*Rm[j,1]+gradN[i,1]*Rm[j,2]
        end
    end
    return B::FFltMat;
end
export Blmat!

function Blmat!(::Type{DeforModelRed2DStress},
                B::FFltMat,N::FFltMat,gradN::FFltMat,c::FFltMat,Rm::FFltMat)
    # Compute the strain-displacement matrix for a two-manifold element.
    #
    # B = Blmat2(self,N,gradN,c,Rm)
    #
    # Compute the linear, displacement independent, strain-displacement matrix
    # for a two-manifold element.
    #
    # Arguments
    #      self=finite element model
    #      N= matrix of basis function values
    #      gradN=matrix of basis function gradients with respect to the
    #         Cartesian coordinates in the directions of the material orientation
    #      c=array of spatial coordinates of the evaluation point
    #         in the global Cartesian coordinates
    #      Rm=orthogonal matrix with the unit basis vectors of the local
    #         material orientation coordinate system as columns;
    #         supply Rm as empty when there is no need for global-to-local
    #         transformation (for instance for isotropic materials).
    #         size(Rm)= [ndim,2], where ndim = number of spatial dimensions
    #         of the embedding space, greater than or equal to  2.
    # Output:
    #      B = strain-displacement matrix, size (B) = [3,nnodes*2],
    #         where  nnodes = number of
    #         finite element nodes on the element.
    #
    nfn= size(gradN,1);
    dim =size(c,2);
    @assert size(B)==(3,dim*nfn)
    fill!(B,0.0);
    @inbounds for i= 1:nfn
        k=dim*(i-1);
        @inbounds for j=1:dim
            B[1,k+j]=gradN[i,1]*Rm[j,1]
            B[2,k+j]=gradN[i,2]*Rm[j,2]
            B[3,k+j]=gradN[i,2]*Rm[j,1]+gradN[i,1]*Rm[j,2]
        end
    end
    return B::FFltMat;
end
export Blmat!

function Blmat!(::Type{DeforModelRed2DAxisymm},
                B::FFltMat,N::FFltMat,gradN::FFltMat,c::FFltMat,Rm::FFltMat)
    # Compute the strain-displacement matrix for a two-manifold element for
    # axially symmetric models.
    #
    #      self=finite element model
    #      N= matrix of basis function values
    #      gradN=matrix of basis function gradients with respect to the
    #         Cartesian coordinates in the directions of the material orientation
    #      c=array of spatial coordinates of the evaluation point
    #         in the global Cartesian coordinates
    #      Rm=orthogonal matrix with the unit basis vectors of the local
    #         material orientation coordinate system as columns;
    #         supply Rm as empty when there is no need for global-to-local
    #         transformation (for instance for isotropic materials).
    #         size(Rm)= [ndim,2], where ndim = number of spatial dimensions
    #         of the embedding space, greater than or equal to  2.
    # Output:
    #      B = strain-displacement matrix, size (B) = [4,nnodes*2],
    #         where  nnodes = number of
    #         finite element nodes on the element.
    #
    nfn= size(gradN,1);
    r=c[1];
    if r <= eps(1.0)
        r = eps(1.0);
    end
    dim =size(c,2);
    @assert size(B)==(4,dim*nfn)
    fill!(B,0.0);
    @inbounds for i= 1:nfn
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
export Blmat!

function Blmat!(::Type{DeforModelRed3D},
                B::FFltMat,N::FFltMat,gradN::FFltMat,c::FFltMat,Rm::FFltMat)
    # Compute the strain-displacement matrix for a three-manifold element.
    #
    # B = Blmat3(self,N,gradN,c,Rm)
    #
    # Compute the linear, displacement independent, strain-displacement matrix
    # for a three-manifold element.
    #
    # Arguments
    #      self=finite element model
    #      N= matrix of basis function values
    #      gradN=matrix of basis function gradients with respect to the
    #         Cartesian coordinates in the directions of the material orientation
    #      c=array of spatial coordinates of the evaluation point
    #         in the global Cartesian coordinates
    #      Rm=orthogonal matrix with the unit basis vectors of the local
    #         material orientation coordinate system as columns;
    #         supply Rm as empty when there is no need for global-to-local
    #         transformation (for instance for isotropic materials).
    #         size(Rm)= [ndim,ndim], where ndim = number of spatial dimensions
    #         of the embedding space.
    # Output:
    #      B = strain-displacement matrix, size (B) = [nstrain,nnodes*dim],
    #         where nstrain= number of strains, dim = Number of spatial
    #         dimensions of the embedding space, and nnodes = number of
    #         finite element nodes on the element.
    #
    nfn= size(gradN,1);
    dim =size(c,2);
    @assert size(B)==(6,dim*nfn)
    fill!(B,0.0);
    @inbounds for i= 1:nfn
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
export Blmat!


end
