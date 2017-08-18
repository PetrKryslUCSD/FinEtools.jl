"""
    FESetModule

Module for finite element sets.
"""
module FESetModule

import Base.count
import Base.cat

export FESet,  FESet0Manifold,  FESet1Manifold,  FESet2Manifold,  FESet3Manifold
export manifdim, nodesperelem, count, getconn!, setlabel!, subset, cat, updateconn!
export FESetP1
export FESetL2, FESetL3
export FESetT3, FESetQ4, FESetQ9, FESetQ8, FESetT6
export FESetH8, FESetH20, FESetH27, FESetT4, FESetT10

using FinEtools.FTypesModule
using FinEtools.RotationUtilModule

abstract type FESet end
abstract type FESet0Manifold <: FESet end
abstract type FESet1Manifold <: FESet end
abstract type FESet2Manifold <: FESet end
abstract type FESet3Manifold <: FESet end

macro add_FESet_fields()
  return esc(:(
  nodesperelem::FInt;
  conn::FIntMat;
  label::FIntVec; )
  )
end
# show(macroexpand(:(@add_FESet_fields)))

"""
    manifdim(me)
Get the manifold dimension.
"""
manifdim(me::FESet0Manifold)::FInt  =0
manifdim(me::FESet1Manifold)::FInt  =1
manifdim(me::FESet2Manifold)::FInt  =2
manifdim(me::FESet3Manifold)::FInt  =3

"""
    nodesperelem(self::T)::FInt where {T<:FESet}

Get the number of nodes  connected  by  the finite element.
"""
(nodesperelem(self::T)::FInt) where {T<:FESet} = self.nodesperelem::FInt

"""
    count(self::T)::FInt where {T<:FESet}

Get the number of individual connectivities in the FE set.
"""
(count(self::T)::FInt) where {T<:FESet} = size(self.conn, 1)

"""
    getconn!(self::T, conn::CC, j::FInt) where {T<:FESet, CC}

Get the connectivity of the jth element in the set.

The connectivity conn[j, :] lists the node  numbers  of the nodes connected
by the jth element.
"""
function getconn!(self::T, conn::CC, j::FInt) where {T<:FESet, CC}
  for i=1:size(self.conn, 2)
    conn[i]=self.conn[j, i];
  end
  return self
end

"""
    boundaryconn(self::T) where {T<:FESet}

Get boundary connectivity.
"""
function boundaryconn(self::T) where {T<:FESet}
    return privboundaryconn(self)
end

"""
    boundaryfe(self::T) where {T<:FESet}

Return the constructor of the type of the boundary finite element.
"""
function boundaryfe(self::T) where {T<:FESet}
    return privboundaryfe(self)
end

"""
    bfun(self::T,  param_coords::FFltVec)::FFltMat where {T<:FESet}

Compute the values of the basis functions at a given parametric coordinate.
"""
function bfun(self::T,  param_coords::FFltVec)::FFltMat where {T<:FESet}
    return privbfun(self, param_coords)
end

"""
    bfundpar(self::T,  param_coords::FFltVec)::FFltMat where {T<:FESet}

Compute the values of the basis function gradients at a given parametric coordinate.
"""
function bfundpar(self::T,  param_coords::FFltVec)::FFltMat where {T<:FESet}
  return privbfundpar(self, param_coords)
end

"""
    Jacobianmatrix!(self::T,  J::FFltMat,  x::FFltMat,
       gradNparams::FFltMat) where {T<:FESet}

Evaluate the Jacobian matrix.
"""
function Jacobianmatrix!(self::T,  J::FFltMat,  x::FFltMat,
   gradNparams::FFltMat) where {T<:FESet}
  @assert (size(J, 2) == size(x, 2)) && (size(J, 1) == size(gradNparams, 2))
  @inbounds for nx = 1:size(x, 2)
    @inbounds for mx = 1:size(gradNparams, 2)
      accumulator::FFlt = 0.0;
      @inbounds for kx = 1:size(x, 1)
        accumulator = accumulator + x[kx, mx]*gradNparams[kx, nx]; # = x' * gradNparams
      end
      J[mx, nx] = accumulator
    end
  end
  return me
end

"""
    setlabel!(self::T, val::FInt) where {T<:FESet}

Set the label of the entire finite elements set.
"""
function setlabel!(self::T, val::FInt) where {T<:FESet}
  self.label = zeros(FInt, size(self.conn, 1));
  fill!(self.label, val)
end

"""
    setlabel!(self::T, val::FIntVec) where {T<:FESet}

Set the labels of individual elements.
"""
function setlabel!(self::T, val::FIntVec) where {T<:FESet}
  #    Set the label of this set.
  @assert size(self.conn, 1)==length(val) "Must get one  label per finite element connectivity"
  self.label=zeros(FInt, size(self.conn, 1));
  copy!(self.label, val);
end

"""
    subset(self::T, L::FIntVec) where {T<:FESet}

Extract a subset of the finite elements from the given finite element set.
"""
function subset(self::T, L::FIntVec) where {T<:FESet}
  result = deepcopy(self)
  result.conn = deepcopy(self.conn[L, :])
  result.label = deepcopy(self.label[L])
  return  result
end

"""
    cat(self::T,  other::T) where {T<:FESet}

Concatenate the connectivities of two FE sets.
"""
function cat(self::T,  other::T) where {T<:FESet}
  @assert self.nodesperelem==other.nodesperelem
  result =deepcopy(self)
  result.conn =[self.conn; other.conn];
  setlabel!(result, vcat(self.label, other.label))
  return result
end

"""
    updateconn!(self::T, NewIDs::FIntVec) where {T<:FESet}

Update the connectivity after the IDs of nodes changed.

NewIDs= new node IDs. Note that indexes in the conn array point
_into_ the  NewIDs array. After the connectivity was updated
this is no longer true!
"""
function updateconn!(self::T, NewIDs::FIntVec) where {T<:FESet}
  for i=1:size(self.conn, 1)
    for j=1:size(self.conn, 2)
      self.conn[i, j]=NewIDs[self.conn[i, j]];
    end
  end
  return self
end

"""
    Jacobian(self::T, J::FFltMat)::FFlt where {T<:FESet1Manifold}

Evaluate the curve Jacobian.

`J` = Jacobian matrix, columns are tangent to parametric coordinates curves.
"""
function Jacobian(self::T, J::FFltMat)::FFlt where {T<:FESet0Manifold}
    return 1.0::FFlt;
end

"""
    Jacobian(self::T, J::FFltMat)::FFlt where {T<:FESet1Manifold}

Evaluate the curve Jacobian.

`J` = Jacobian matrix, columns are tangent to parametric coordinates curves.
"""
function Jacobian(self::T, J::FFltMat)::FFlt where {T<:FESet1Manifold}
  sdim,  ntan = size(J);
  @assert ntan == 1 "Expected number of tangent vectors: 1"
  return norm(J)::FFlt;
end

"""
    gradN!(self::FESet1Manifold, gradN::FFltMat, gradNparams::FFltMat,
      redJ::FFltMat)

Compute the gradient of the basis functions with the respect to
the "reduced" spatial coordinates.

gradN= output,  matrix of gradients,  one per row
gradNparams= matrix of gradients with respect to parametric coordinates, one per row
redJ= reduced Jacobian matrix `redJ=transpose(Rm)*J`
"""
function gradN!(self::FESet1Manifold, gradN::FFltMat, gradNparams::FFltMat, redJ::FFltMat)
  for r=1:size(gradN, 1)
    gradN[r, 1]=  gradNparams[r, 1]/redJ[1, 1]
  end
end

"""
    Jacobian(self::T, J::FFltMat)::FFlt where {T<:FESet2Manifold}

Evaluate the curve Jacobian.

`J` = Jacobian matrix, columns are tangent to parametric coordinates curves.
"""
function Jacobian(self::T, J::FFltMat)::FFlt where {T<:FESet2Manifold}
  sdim,  ntan = size(J);
  @assert ntan == 2 "Expected number of tangent vectors: 2"
  if sdim == ntan
    @inbounds Jac = (J[1, 1]*J[2, 2] - J[2, 1]*J[1, 2])
    return Jac::FFlt;# is det(J);% Compute the Jacobian
  else
    return norm(RotationUtilModule.cross(J[:, 1], J[:, 2]))::FFlt;
  end
end

"""
    gradN!(self::FESet2Manifold, gradN::FFltMat, gradNparams::FFltMat,
      redJ::FFltMat)

Compute the gradient of the basis functions with the respect to
the "reduced" spatial coordinates.

gradN= output,  matrix of gradients,  one per row
gradNparams= matrix of gradients with respect to parametric coordinates, one per row
redJ= reduced Jacobian matrix `redJ=transpose(Rm)*J`
"""
function gradN!(self::FESet2Manifold, gradN::FFltMat, gradNparams::FFltMat, redJ::FFltMat)
  # This is the unrolled version that avoids allocation of a 3 x 3 matrix
  invdet=1.0/(redJ[1, 1]*redJ[2, 2] - redJ[1, 2]*redJ[2, 1]);
  invredJ11 =  (redJ[2, 2])*invdet;
  invredJ12 = -(redJ[1, 2])*invdet;
  invredJ21 = -(redJ[2, 1])*invdet;
  invredJ22 =  (redJ[1, 1])*invdet;
  @assert size(gradN, 1)==size(gradNparams, 1)
  @inbounds for r=1:size(gradN, 1)
    gradN[r, 1]= gradNparams[r, 1]*invredJ11 +gradNparams[r, 2]*invredJ21;
    gradN[r, 2]= gradNparams[r, 1]*invredJ12 +gradNparams[r, 2]*invredJ22;
  end
end

"""
    Jacobian(self::T, J::FFltMat)::FFlt where {T<:FESet3Manifold}

Evaluate the volume Jacobian.

`J` = Jacobian matrix, columns are tangent to parametric coordinates curves.
"""
function Jacobian(self::T, J::FFltMat)::FFlt where {T<:FESet3Manifold}
  sdim,  ntan = size(J);
  @assert (ntan == 3) && (sdim == 3) "Expected number of tangent vectors: 3"
  #Jac = det(J);# Compute the Jacobian
  # The unrolled version
  return (+J[1, 1]*(J[2, 2]*J[3, 3]-J[3, 2]*J[2, 3])
          -J[1, 2]*(J[2, 1]*J[3, 3]-J[2, 3]*J[3, 1])
          +J[1, 3]*(J[2, 1]*J[3, 2]-J[2, 2]*J[3, 1]) )::FFlt;
end

"""
    gradN!(self::FESet3Manifold, gradN::FFltMat, gradNparams::FFltMat,
      redJ::FFltMat)

Compute the gradient of the basis functions with the respect to
the "reduced" spatial coordinates.

gradN= output,  matrix of gradients,  one per row
gradNparams= matrix of gradients with respect to parametric coordinates, one per row
redJ= reduced Jacobian matrix `redJ=transpose(Rm)*J`
"""
function gradN!(self::FESet3Manifold, gradN::FFltMat, gradNparams::FFltMat, redJ::FFltMat)
  invdet = 1.0 / ( +redJ[1, 1]*(redJ[2, 2]*redJ[3, 3]-redJ[3, 2]*redJ[2, 3])
  -redJ[1, 2]*(redJ[2, 1]*redJ[3, 3]-redJ[2, 3]*redJ[3, 1])
  +redJ[1, 3]*(redJ[2, 1]*redJ[3, 2]-redJ[2, 2]*redJ[3, 1]) );
  # This is the unrolled version that avoids allocation of a 3 x 3 matrix
  invredJ11 =  (redJ[2, 2]*redJ[3, 3]-redJ[3, 2]*redJ[2, 3])*invdet;
  invredJ12 = -(redJ[1, 2]*redJ[3, 3]-redJ[1, 3]*redJ[3, 2])*invdet;
  invredJ13 =  (redJ[1, 2]*redJ[2, 3]-redJ[1, 3]*redJ[2, 2])*invdet;
  invredJ21 = -(redJ[2, 1]*redJ[3, 3]-redJ[2, 3]*redJ[3, 1])*invdet;
  invredJ22 =  (redJ[1, 1]*redJ[3, 3]-redJ[1, 3]*redJ[3, 1])*invdet;
  invredJ23 = -(redJ[1, 1]*redJ[2, 3]-redJ[2, 1]*redJ[1, 3])*invdet;
  invredJ31 =  (redJ[2, 1]*redJ[3, 2]-redJ[3, 1]*redJ[2, 2])*invdet;
  invredJ32 = -(redJ[1, 1]*redJ[3, 2]-redJ[3, 1]*redJ[1, 2])*invdet;
  invredJ33 =  (redJ[1, 1]*redJ[2, 2]-redJ[2, 1]*redJ[1, 2])*invdet;
  @assert size(gradN, 1)==size(gradNparams, 1)
  @inbounds for r=1:size(gradN, 1)
    gradN[r, 1]= gradNparams[r, 1]*invredJ11 +gradNparams[r, 2]*invredJ21 +gradNparams[r, 3]*invredJ31;
    gradN[r, 2]= gradNparams[r, 1]*invredJ12 +gradNparams[r, 2]*invredJ22 +gradNparams[r, 3]*invredJ32;
    gradN[r, 3]= gradNparams[r, 1]*invredJ13 +gradNparams[r, 2]*invredJ23 +gradNparams[r, 3]*invredJ33;
  end
end



"""
    FESetP1

Type for sets of point-like of finite elements.
"""
mutable struct FESetP1 <: FESet0Manifold
  @add_FESet_fields

  function    FESetP1(conn::FIntMat=[])
    nodesperelem::FInt  = 1
    @assert (size(conn, 2) == nodesperelem) "Number of nodes per element mismatched"
    # Need to make a COPY of the input arrays
    self =new(nodesperelem, deepcopy(conn), deepcopy(FInt[]))
    setlabel!(self, 0)
    return self
  end
end

privbfun(self::FESetP1,  param_coords::FFltVec) = reshape([1.0], 1, 1) # make sure this is a matrix
privbfundpar(self::FESetP1,  param_coords::FFltVec) = zeros(1, 0)

function privboundaryconn(self::FESetP1)
    # Get boundary connectivity.
    return [];
end

function privboundaryfe(self::FESetP1)
    return None;
end

"""
    FESetL2

Type for sets of curve-like of finite elements with two nodes.
"""
mutable struct FESetL2 <: FESet1Manifold
  @add_FESet_fields

  function    FESetL2(conn::FIntMat=[])
    nodesperelem::FInt  = 2
    @assert (size(conn, 2) == nodesperelem) "Number of nodes per element mismatched"
    # Need to make a COPY of the input arrays
    self =new(nodesperelem, deepcopy(conn), deepcopy(FInt[]))
    setlabel!(self, 0)
    return self
  end
end

privbfun(self::FESetL2,  param_coords::FFltVec) = reshape([(1. - param_coords[1]); (1. + param_coords[1])] / 2.0, 2, 1) # make sure this is a matrix
privbfundpar(self::FESetL2,  param_coords::FFltVec) = reshape([-1.0; +1.0]/2.0, 2, 1)

function privboundaryconn(self::FESetL2)
    # Get boundary connectivity.
    return [self.conn[:, 1];self.conn[:, 2]];
end

function privboundaryfe(self::FESetL2)
    return FESetP1;
end

"""
    FESetL3

Type for sets of curve-like of finite elements with three nodes.
"""
mutable struct FESetL3 <: FESet1Manifold
  @add_FESet_fields

  function    FESetL3(conn::FIntMat=[])
    nodesperelem::FInt  = 3
    @assert (size(conn, 2) == nodesperelem) "Number of nodes per element mismatched"
    # Need to make a COPY of the input arrays
    self =new(nodesperelem, deepcopy(conn), deepcopy(FInt[]))
    setlabel!(self, 0)
    return self
  end
end

function privbfun(self::FESetL3,  param_coords::FFltVec)
    xi=param_coords[1];
    val = [(xi-1)*xi/2;  (xi+1)*xi/2;  -(xi+1)*(xi-1)];
    return reshape(val, 3, 1) # make sure this is a matrix
end

function privbfundpar(self::FESetL3,  param_coords::FFltVec)
    xi=param_coords[1];
    val = [(xi-1/2); (xi+1/2); -2*xi];
    return reshape(val, 3, 1) # make sure this is a matrix
end

function privboundaryconn(self::FESetL3)
    return [self.conn[:, 1];self.conn[:, 2]];
end

function privboundaryfe(self::FESetL3)
    return FESetP1;
end

"""
    FESetT3

Type for sets of surface-like of triangular finite elements with three nodes.
"""
mutable struct FESetT3 <: FESet2Manifold
  @add_FESet_fields

  function  FESetT3( conn::FIntMat=[])
    nodesperelem::FInt  = 3
    @assert (size(conn, 2) == nodesperelem) "Number of nodes per element mismatched"
    # Need to make a COPY of the input arrays
    self = new(nodesperelem, deepcopy(conn), deepcopy(FInt[]))
    setlabel!(self, 0)
    return self
  end
end

function privbfun(self::FESetT3,  param_coords::FFltVec)
    # Evaluate the basis function matrix for an 3-node triangle.
    return reshape([(1 - param_coords[1] - param_coords[2]); param_coords[1]; param_coords[2]], 3, 1); # Make sure this is a matrix
end

function privbfundpar(self::FESetT3,  param_coords::FFltVec)
    # Evaluate the derivatives of the basis function matrix.
    return  [-1. -1.;  +1.  0.;  0. +1.];
end

function privboundaryconn(self::FESetT3)
    return [self.conn[:, 1:2];self.conn[:, 2:3];self.conn[:, [3, 1]]];
end

function privboundaryfe(self::FESetT3)
    # Get  the constructor of the class of the  boundary finite element.
    return FESetL2;
end


"""
    FESetQ4

Type for sets of surface-like of quadrilateral finite elements with four nodes.
"""
mutable struct FESetQ4 <: FESet2Manifold
  @add_FESet_fields

  function    FESetQ4(conn::FIntMat=[])
    nodesperelem::FInt  = 4
    @assert (size(conn, 2) == nodesperelem) "Number of nodes per element mismatched"
    # Need to make a COPY of the input arrays
    self =new(nodesperelem, deepcopy(conn), deepcopy(FInt[]))
    setlabel!(self, 0)
    return self
  end
end

function privbfun(self::FESetQ4,  param_coords::FFltVec)
    # Evaluate the basis function matrix for an 4-node quadrilateral.
    val = [0.25 * (1. - param_coords[1]) * (1. - param_coords[2]);
           0.25 * (1. + param_coords[1]) * (1. - param_coords[2]);
           0.25 * (1. + param_coords[1]) * (1. + param_coords[2]);
           0.25 * (1. - param_coords[1]) * (1. + param_coords[2])];
    return reshape(val, 4, 1); # Make sure this is a matrix
end

function privbfundpar(self::FESetQ4,  param_coords::FFltVec)
    # Evaluate the derivatives of the basis function matrix.
    val = [-(1. - param_coords[2])*0.25 -(1. - param_coords[1])*0.25;
            (1. - param_coords[2])*0.25 -(1. + param_coords[1])*0.25;
            (1. + param_coords[2])*0.25 (1. + param_coords[1])*0.25;
           -(1. + param_coords[2])*0.25 (1. - param_coords[1])*0.25];
    return val
end

function privboundaryconn(self::FESetQ4)
    return [self.conn[:, 1:2];self.conn[:, 2:3];self.conn[:, 3:4];self.conn[:, [4, 1]]];
end

function privboundaryfe(self::FESetQ4)
    # Get  the constructor of the class of the  boundary finite element.
    return FESetL2;
end

"""
    FESetQ9

Type for sets of surface-like of quadrilateral finite elements with nine nodes.
"""
mutable struct FESetQ9 <: FESet2Manifold
  @add_FESet_fields

  function    FESetQ9(conn::FIntMat=[])
    nodesperelem::FInt  = 9
    @assert (size(conn, 2) == nodesperelem) "Number of nodes per element mismatched"
    # Need to make a COPY of the input arrays
    self =new(nodesperelem, deepcopy(conn), deepcopy(FInt[]))
    setlabel!(self, 0)
    return self
  end
end

function privbfun(self::FESetQ9,  param_coords::FFltVec)
    # Evaluate the basis function matrix for an 4-node quadrilateral.
    xi=param_coords[1];
    xis = [(xi-1)*xi/2;  (xi+1)*xi/2;  -(xi+1)*(xi-1)];
    eta=param_coords[2];
    etas = [(eta-1)*eta/2;  (eta+1)*eta/2;  -(eta+1)*(eta-1)];
    xisetas=(xis*etas');
    val = xisetas([1     2     5     4     3     8     6     7     9])';
    return reshape(val, 9, 1); # Make sure this is a matrix
end

function privbfundpar(self::FESetQ9,  param_coords::FFltVec)
    # Evaluate the derivatives of the basis function matrix.
    xi=param_coords[1];
    xis = [(xi-1)*xi/2;  (xi+1)*xi/2;  -(xi+1)*(xi-1)];
    dxis = [(xi-1/2); (xi+1/2); -2*xi];
    eta=param_coords[2];
    etas = [(eta-1)*eta/2;  (eta+1)*eta/2;  -(eta+1)*(eta-1)];
    detas=[(eta-1/2); (eta+1/2); -2*eta];
    dxisetas=(dxis*etas');
    xisdetas=(xis*detas');
    ix=[1     2     5     4     3     8     6     7     9]
    val =  [dxisetas(ix)'             xisdetas(ix)'];
    return reshape(val, 9, 2); # Make sure this is a matrix
end

function privboundaryconn(self::FESetQ9)
    return [self.conn[:, [1, 2, 5]];self.conn[:, [2, 3, 6]];self.conn[:, [3, 4, 7]];self.conn[:, [4, 1, 8]]];
end

function privboundaryfe(self::FESetQ9)
    return FESetL3;
end

"""
    FESetQ8

Type for sets of surface-like of quadrilateral finite elements with eight nodes.
"""
mutable struct FESetQ8 <: FESet2Manifold
  @add_FESet_fields

  function    FESetQ8(conn::FIntMat=[])
    nodesperelem::FInt  = 8
    @assert (size(conn, 2) == nodesperelem) "Number of nodes per element mismatched"
    # Need to make a COPY of the input arrays
    self =new(nodesperelem, deepcopy(conn), deepcopy(FInt[]))
    setlabel!(self, 0)
    return self
  end
end

function privbfun(self::FESetQ8,  param_coords::FFltVec)
    # Evaluate the basis function matrix for an 4-node quadrilateral.
    xim    = (-1 + param_coords[1]);
    etam   = (-1 + param_coords[2]);
    xip    = (1 + param_coords[1]);
    etap   = (1 + param_coords[2]);
    val = [ -1.0/4*xim*etam*(1+param_coords[1]+param_coords[2]);
           1.0/4*xip*etam*(1-param_coords[1]+param_coords[2]);
           -1.0/4*xip*etap*(1-param_coords[1]-param_coords[2]);
           1.0/4*xim*etap*(1+param_coords[1]-param_coords[2]);
           1.0/2*xim*xip*etam;
           -1.0/2*etam*etap*xip;
           -1.0/2*xim*xip*etap;
           1.0/2*etam*etap*xim];
    return reshape(val, 8, 1); # Make sure this is a matrix
end

function privbfundpar(self::FESetQ8,  param_coords::FFltVec)
    # Evaluate the derivatives of the basis function matrix.
    xi =param_coords[1];
    eta =param_coords[2];
    val=zeros(FFlt, 8, 2)
    val[:, 1]= [1.0/4*(1-eta)*(1+xi+eta)-1.0/4*(1-xi)*(1-eta);
               -1.0/4*(1-eta)*(1-xi+eta)+1.0/4*(1+xi)*(1-eta);
               -1.0/4*(1+eta)*(1-xi-eta)+1.0/4*(1+xi)*(1+eta);
               1.0/4*(1+eta)*(1+xi-eta)-1.0/4*(1-xi)*(1+eta);
               -1.0/2*(1+xi)*(1-eta)+1.0/2*(1-xi)*(1-eta);
               1.0/2*(1-eta)*(1+eta);
               -1.0/2*(1+xi)*(1+eta)+1.0/2*(1-xi)*(1+eta);
               -1.0/2*(1-eta)*(1+eta)];
    val[:, 2] = [1.0/4*(1-xi)*(1+xi+eta)-1.0/4*(1-xi)*(1-eta);
                1.0/4*(1+xi)*(1-xi+eta)-1.0/4*(1+xi)*(1-eta);
                -1.0/4*(1+xi)*(1-xi-eta)+1.0/4*(1+xi)*(1+eta);
                -1.0/4*(1-xi)*(1+xi-eta)+1.0/4*(1-xi)*(1+eta);
                -1.0/2*(1-xi)*(1+xi);
                -1.0/2*(1+xi)*(1+eta)+1.0/2*(1+xi)*(1-eta);
                1.0/2*(1-xi)*(1+xi);
                -1.0/2*(1-xi)*(1+eta)+1.0/2*(1-xi)*(1-eta)];
    return  reshape(val, 8, 2); # Make sure this is a matrix
end

function privboundaryconn(self::FESetQ8)
    # Get boundary connectivity.
    return [self.conn[:, [1, 2, 5]];self.conn[:, [2, 3, 6]];self.conn[:, [3, 4, 7]];self.conn[:, [4, 1, 8]]];
end

function privboundaryfe(self::FESetQ8)
    return FESetL3;
end

"""
    FESetT6

Type for sets of surface-like of triangular finite elements with six nodes.
"""
mutable struct FESetT6 <: FESet2Manifold
  @add_FESet_fields

  function    FESetT6(conn::FIntMat=[])
    nodesperelem::FInt  = 6
    @assert (size(conn, 2) == nodesperelem) "Number of nodes per element mismatched"
    # Need to make a COPY of the input arrays
    self =new(nodesperelem, deepcopy(conn), deepcopy(FInt[]))
    setlabel!(self, 0)
    return self
  end
end

function privbfun(self::FESetT6,  param_coords::FFltVec)
    # Evaluate the basis function matrix for an 4-node quadrilateral.
    r=param_coords[1];
    s=param_coords[2];
    t = 1. - r - s;
    val = [t * (t + t - 1);
           r * (r + r - 1);
           s * (s + s - 1);
           4 * r * t;
           4 * r * s;
           4 * s * t];
    return reshape(val, 6, 1); # Make sure this is a matrix
end

function privbfundpar(self::FESetT6,  param_coords::FFltVec)
    # Evaluate the derivatives of the basis function matrix.
    r =param_coords[1];
    s =param_coords[2];
    t = 1. - r - s;
    val = [-3+4*r+4*s  -3+4*r+4*s;
           4*r-1  0.0;
           0.0 4*s-1;
           4-8*r-4*s  -4*r;
           4*s  4*r;
           -4*s  4-4*r-8*s];
    return  reshape(val, 6, 2); # Make sure this is a matrix
end

function privboundaryconn(self::FESetT6)
    # Get boundary connectivity.
    return [self.conn[:, [1,  2,  4]];self.conn[:, [2,  3,  5]];self.conn[:, [3,  1,  6]]];
end

function privboundaryfe(self::FESetT6)
    # Get  the constructor of the class of the  boundary finite element.
    return FESetL3;
end

"""
    FESetH8

Type for sets of volume-like of hexahedral finite elements with eight nodes.
"""
mutable struct FESetH8 <: FESet3Manifold
  @add_FESet_fields

  function    FESetH8(conn::FIntMat=[])
    nodesperelem::FInt  = 8
    @assert (size(conn, 2) == nodesperelem) "Number of nodes per element mismatched"
    # Need to make a COPY of the input arrays
    self =new(nodesperelem, deepcopy(conn), deepcopy(FInt[]))
    setlabel!(self, 0)
    return self
  end
end

function privbfun(self::FESetH8,  param_coords::FFltVec)
    # Evaluate the basis function matrix for an 8-node hexahedron.
    one_minus_xi    = (1.0 - param_coords[1]);
    one_minus_eta   = (1.0 - param_coords[2]);
    one_minus_theta = (1.0 - param_coords[3]);
    one_plus_xi     = (1.0 + param_coords[1]);
    one_plus_eta    = (1.0 + param_coords[2]);
    one_plus_theta  = (1.0 + param_coords[3]);
    val = [one_minus_xi*one_minus_eta*one_minus_theta;
           one_plus_xi*one_minus_eta*one_minus_theta;
           one_plus_xi*one_plus_eta*one_minus_theta;
           one_minus_xi*one_plus_eta*one_minus_theta;
           one_minus_xi*one_minus_eta*one_plus_theta;
           one_plus_xi*one_minus_eta*one_plus_theta;
           one_plus_xi*one_plus_eta*one_plus_theta;
           one_minus_xi*one_plus_eta*one_plus_theta] / 8.0;
    return reshape(val, 8, 1); # Make sure this is a matrix
end

function privbfundpar(self::FESetH8,  param_coords::FFltVec)
    # Evaluate the derivatives of the basis function matrix.
     omxi    = (1.0 - param_coords[1]);
    ometa   = (1.0 - param_coords[2]);
    omtheta = (1.0 - param_coords[3]);
    opxi     = (1.0 + param_coords[1]);
    opeta    = (1.0 + param_coords[2]);
    optheta  = (1.0 + param_coords[3]);
   val = [-ometa*omtheta  ometa*omtheta  opeta*omtheta  -opeta*omtheta  -ometa*optheta  ometa*optheta  opeta*optheta   -opeta*optheta;
          -omxi*omtheta  -opxi*omtheta   opxi*omtheta   omxi*omtheta    -omxi*optheta   -opxi*optheta  opxi*optheta   omxi*optheta;
          -omxi*ometa    -opxi*ometa     -opxi*opeta    -omxi*opeta      omxi*ometa     opxi*ometa     opxi*opeta    omxi*opeta]'/8.0;
    return val
end

function privboundaryconn(self::FESetH8)
    return  [self.conn[:, [1,  4,  3,  2]];
             self.conn[:, [1,  2,  6,  5]];
             self.conn[:, [2,  3,  7,  6]];
             self.conn[:, [3,  4,  8,  7]];
             self.conn[:, [4,  1,  5,  8]];
             self.conn[:, [6,  7,  8,  5]]];
end

function privboundaryfe(self::FESetH8)
    # Get  the constructor of the class of the  boundary finite element.
    return FESetQ4;
end

"""
    FESetH20

Type for sets of volume-like of hexahedral finite elements with 20 nodes.
"""
mutable struct FESetH20 <: FESet3Manifold
  @add_FESet_fields

  function    FESetH20(conn::FIntMat=[])
    nodesperelem::FInt  = 20
    @assert (size(conn, 2) == nodesperelem) "Number of nodes per element mismatched"
    self =new(nodesperelem, deepcopy(conn), deepcopy(FInt[]))
    setlabel!(self, 0)
    return self
  end
end

function privbfun(self::FESetH20,  param_coords::FFltVec)
    # Evaluate the basis function matrix for an 20-node hexahedron.
    xim    = (-1 + param_coords[1]);
    etam   = (-1 + param_coords[2]);
    zetam  = (-1 + param_coords[3]);
    xip    = (1 + param_coords[1]);
    etap   = (1 + param_coords[2]);
    zetap  = (1 + param_coords[3]);
    val = [ 1.0/8*xim*etam*zetam*(2+param_coords[1]+param_coords[2]+param_coords[3]);
           -1.0/8*xip*etam*zetam*(2-param_coords[1]+param_coords[2]+param_coords[3]);
           1.0/8*xip*etap*zetam*(2-param_coords[1]-param_coords[2]+param_coords[3]);
           -1.0/8*xim*etap*zetam*(2+param_coords[1]-param_coords[2]+param_coords[3]);
           1.0/8*xim*etam*zetap*(-2-param_coords[1]-param_coords[2]+param_coords[3]);
           -1.0/8*xip*etam*zetap*(-2+param_coords[1]-param_coords[2]+param_coords[3]);
           1.0/8*xip*etap*zetap*(-2+param_coords[1]+param_coords[2]+param_coords[3]);
           -1.0/8*xim*etap*zetap*(-2-param_coords[1]+param_coords[2]+param_coords[3]);
           -1.0/4*xim*xip*etam*zetam;
           1.0/4*etam*etap*xip*zetam;
           1.0/4*xim*xip*etap*zetam;
           -1.0/4*etam*etap*xim*zetam;
           1.0/4*xim*xip*etam*zetap;
           -1.0/4*etam*etap*xip*zetap;
           -1.0/4*xim*xip*etap*zetap;
           1.0/4*etam*etap*xim*zetap;
           -1.0/4*zetam*zetap*xim*etam;
           1.0/4*zetam*zetap*xip*etam;
           -1.0/4*zetam*zetap*xip*etap;
           1.0/4*zetam*zetap*xim*etap];
    return reshape(val, 20, 1); # Make sure this is a matrix
end

function privbfundpar(self::FESetH20,  param_coords::FFltVec)
    # Evaluate the derivatives of the basis function matrix.
    xim    = -(-1 + param_coords[1]);
    etam   = -(-1 + param_coords[2]);
    zetam  = -(-1 + param_coords[3]);
    xip    = (1 + param_coords[1]);
    etap   = (1 + param_coords[2]);
    zetap  = (1 + param_coords[3]);
    twoppp =(2+param_coords[1]+param_coords[2]+param_coords[3]);
    twompp =(2-param_coords[1]+param_coords[2]+param_coords[3]);
    twopmp =(2+param_coords[1]-param_coords[2]+param_coords[3]);
    twoppm =(2+param_coords[1]+param_coords[2]-param_coords[3]);
    twommp =(2-param_coords[1]-param_coords[2]+param_coords[3]);
    twopmm =(2+param_coords[1]-param_coords[2]-param_coords[3]);
    twompm =(2-param_coords[1]+param_coords[2]-param_coords[3]);
    twommm =(2-param_coords[1]-param_coords[2]-param_coords[3]);
    val =zeros(FFlt, 20, 3)
   val[:, 1]= [1.0/8*etam*zetam*twoppp-1.0/8*xim*etam*zetam;
           -1.0/8*etam*zetam*twompp+1.0/8*xip*etam*zetam;
           -1.0/8*etap*zetam*twommp+1.0/8*xip*etap*zetam;
           1.0/8*etap*zetam*twopmp-1.0/8*xim*etap*zetam;
           1.0/8*etam*zetap*twoppm-1.0/8*xim*etam*zetap;
           -1.0/8*etam*zetap*twompm+1.0/8*xip*etam*zetap;
           -1.0/8*etap*zetap*twommm+1.0/8*xip*etap*zetap;
           1.0/8*etap*zetap*twopmm-1.0/8*xim*etap*zetap;
           -1.0/4*xip*etam*zetam+1.0/4*xim*etam*zetam;
           1.0/4*etam*etap*zetam;
           -1.0/4*xip*etap*zetam+1.0/4*xim*etap*zetam;
           -1.0/4*etam*etap*zetam;
           -1.0/4*xip*etam*zetap+1.0/4*xim*etam*zetap;
           1.0/4*etam*etap*zetap;
           -1.0/4*xip*etap*zetap+1.0/4*xim*etap*zetap;
           -1.0/4*etam*etap*zetap;
           -1.0/4*zetam*zetap*etam;
           1.0/4*zetam*zetap*etam;
           1.0/4*zetam*zetap*etap;
           -1.0/4*zetam*zetap*etap];
    val[:, 2]= [1.0/8*xim*zetam*twoppp-1.0/8*xim*etam*zetam;
             1.0/8*xip*zetam*twompp-1.0/8*xip*etam*zetam;
             -1.0/8*xip*zetam*twommp+1.0/8*xip*etap*zetam;
             -1.0/8*xim*zetam*twopmp+1.0/8*xim*etap*zetam;
             1.0/8*xim*zetap*twoppm-1.0/8*xim*etam*zetap;
             1.0/8*xip*zetap*twompm-1.0/8*xip*etam*zetap;
             -1.0/8*xip*zetap*twommm+1.0/8*xip*etap*zetap;
             -1.0/8*xim*zetap*twopmm+1.0/8*xim*etap*zetap;
             -1.0/4*xim*xip*zetam;
             -1.0/4*xip*etap*zetam+1.0/4*xip*etam*zetam;
             1.0/4*xim*xip*zetam;
             -1.0/4*xim*etap*zetam+1.0/4*xim*etam*zetam;
             -1.0/4*xim*xip*zetap;
             -1.0/4*xip*etap*zetap+1.0/4*xip*etam*zetap;
             1.0/4*xim*xip*zetap;
             -1.0/4*xim*etap*zetap+1.0/4*xim*etam*zetap;
             -1.0/4*zetam*zetap*xim;
             -1.0/4*zetam*zetap*xip;
             1.0/4*zetam*zetap*xip;
             1.0/4*zetam*zetap*xim];
    val[:, 3]= [1.0/8*xim*etam*twoppp-1.0/8*xim*etam*zetam;
             1.0/8*xip*etam*twompp-1.0/8*xip*etam*zetam;
             1.0/8*xip*etap*twommp-1.0/8*xip*etap*zetam;
             1.0/8*xim*etap*twopmp-1.0/8*xim*etap*zetam;
             -1.0/8*xim*etam*twoppm+1.0/8*xim*etam*zetap;
             -1.0/8*xip*etam*twompm+1.0/8*xip*etam*zetap;
             -1.0/8*xip*etap*twommm+1.0/8*xip*etap*zetap;
             -1.0/8*xim*etap*twopmm+1.0/8*xim*etap*zetap;
             -1.0/4*xim*xip*etam;
             -1.0/4*etam*etap*xip;
             -1.0/4*xim*xip*etap;
             -1.0/4*etam*etap*xim;
             1.0/4*xim*xip*etam;
             1.0/4*etam*etap*xip;
             1.0/4*xim*xip*etap;
             1.0/4*etam*etap*xim;
             -1.0/4*xim*etam*zetap+1.0/4*xim*etam*zetam;
             -1.0/4*xip*etam*zetap+1.0/4*xip*etam*zetam;
             -1.0/4*xip*etap*zetap+1.0/4*xip*etap*zetam;
             -1.0/4*xim*etap*zetap+1.0/4*xim*etap*zetam];
    return reshape(val, 20, 3);
end

function privboundaryconn(self::FESetH20)
    return  [self.conn[:, [1,  4,  3,  2,  12,  11,  10,  9]];
             self.conn[:, [1,  2,  6,  5,  9,  18,  13,  17]];
             self.conn[:, [2,  3,  7,  6,  10,  19,  14,  18]];
             self.conn[:, [3,  4,  8,  7,   11,  20,  15,  19]];
             self.conn[:, [4,  1,  5,  8,  12,  17,  16,  20]];
             self.conn[:, [6,  7,  8,  5,  14,  15,  16,  13]]];
end

function privboundaryfe(self::FESetH20)
    # Get  the constructor of the class of the  boundary finite element.
    return FESetQ8;
end

"""
    FESetH27

Type for sets of volume-like of hexahedral finite elements with 27 nodes.
"""
mutable struct FESetH27 <: FESet3Manifold
  @add_FESet_fields

  function    FESetH27(conn::FIntMat=[])
    nodesperelem::FInt  = 27
    @assert (size(conn, 2) == nodesperelem) "Number of nodes per element mismatched"
    self =new(nodesperelem, deepcopy(conn), deepcopy(FInt[]))
    setlabel!(self, 0)
    return self
  end
end

function privbfun(self::FESetH27,  param_coords::FFltVec)
    # Evaluate the basis function matrix for an 8-node hexahedron.
    xi=param_coords[1];
    eta=param_coords[2];
    zet=param_coords[3];
    x1 =(xi-1);
    y1=(eta-1);
    z1 =(zet-1);
    x2 =(xi+1);
    y2=(eta+1);
    z2 =(zet+1);
    val =[1.0/8.0*z1*zet*x1*xi*y1*eta
          1.0/8.0*z1*zet*x2*xi*y1*eta
          1.0/8.0*z1*zet*x2*xi*y2*eta
          1.0/8.0*z1*zet*x1*xi*y2*eta
          1.0/8.0*z2*zet*x1*xi*y1*eta
          1.0/8.0*z2*zet*x2*xi*y1*eta
          1.0/8.0*z2*zet*x2*xi*y2*eta
          1.0/8.0*z2*zet*x1*xi*y2*eta
          1.0/4.0*z1*zet*(-x2)*x1*y1*eta
          1.0/4.0*z1*zet*x2*xi*(-y2)*y1
          1.0/4.0*z1*zet*(-x2)*x1*y2*eta
          1.0/4.0*z1*zet*x1*xi*(-y2)*y1
          1.0/4.0*z2*zet*(-x2)*x1*y1*eta
          1.0/4.0*z2*zet*x2*xi*(-y2)*y1
          1.0/4.0*z2*zet*(-x2)*x1*y2*eta
          1.0/4.0*z2*zet*x1*xi*(-y2)*y1
          1.0/4.0*(-z2)*z1*x1*xi*y1*eta
          1.0/4.0*(-z2)*z1*x2*xi*y1*eta
          1.0/4.0*(-z2)*z1*x2*xi*y2*eta
          1.0/4.0*(-z2)*z1*x1*xi*y2*eta
          1.0/2.0*z1*zet*(-x2)*x1*(-y2)*y1
          1.0/2.0*(-z2)*z1*(-x2)*x1*y1*eta
          1.0/2.0*(-z2)*z1*x2*xi*(-y2)*y1
          1.0/2.0*(-z2)*z1*(-x2)*x1*y2*eta
          1.0/2.0*(-z2)*z1*x1*xi*(-y2)*y1
          1.0/2.0*z2*zet*(-x2)*x1*(-y2)*y1
          (-z2)*z1*(-x2)*x1*(-y2)*y1];
    return reshape(val, 27, 1); # Make sure this is a matrix
end

function privbfundpar(self::FESetH27,  param_coords::FFltVec)
    # Evaluate the derivatives of the basis function matrix.
    xi=param_coords[1];
    eta=param_coords[2];
    zet=param_coords[3];
    x1 =(xi-1.0/2.0);
    x2 =(xi+1.0/2.0);
    x3 =(xi-1.0);
    x4 =(xi+1.0);
    z1 =(zet-1.0);
    z2 =(zet-1.0/2.0);
    z3 =(zet+1.0);
    z4 =(zet+1.0/2.0);
    y1 =(eta-1.0);
    y2 =(eta-1.0/2.0);
    y3 =(eta+1.0);
    y4 = (eta+1.0/2.0);
    val = [      1.0/4.0*z1*zet*x1*y1*eta        1.0/4.0*z1*zet*x3*xi*y2        1.0/4.0*z2*x3*xi*y1*eta;
                  1.0/4.0*z1*zet*x2*y1*eta        1.0/4.0*z1*zet*x4*xi*y2        1.0/4.0*z2*x4*xi*y1*eta;
                  1.0/4.0*z1*zet*x2*y3*eta        1.0/4.0*z1*zet*x4*xi*y4        1.0/4.0*z2*x4*xi*y3*eta;
                  1.0/4.0*z1*zet*x1*y3*eta        1.0/4.0*z1*zet*x3*xi*y4        1.0/4.0*z2*x3*xi*y3*eta;
                  1.0/4.0*z3*zet*x1*y1*eta        1.0/4.0*z3*zet*x3*xi*y2        1.0/4.0*z4*x3*xi*y1*eta;
                  1.0/4.0*z3*zet*x2*y1*eta        1.0/4.0*z3*zet*x4*xi*y2        1.0/4.0*z4*x4*xi*y1*eta;
                  1.0/4.0*z3*zet*x2*y3*eta        1.0/4.0*z3*zet*x4*xi*y4        1.0/4.0*z4*x4*xi*y3*eta;
                  1.0/4.0*z3*zet*x1*y3*eta        1.0/4.0*z3*zet*x3*xi*y4        1.0/4.0*z4*x3*xi*y3*eta;
                       -1.0/2.0*z1*zet*xi*y1*eta   1.0/2.0*z1*zet*(-x4)*x3*y2   1.0/2.0*z2*(-x4)*x3*y1*eta;
             1.0/2.0*z1*zet*x2*(-y3)*y1             -1.0/2.0*z1*zet*x4*xi*eta   1.0/2.0*z2*x4*xi*(-y3)*y1;
                       -1.0/2.0*z1*zet*xi*y3*eta   1.0/2.0*z1*zet*(-x4)*x3*y4   1.0/2.0*z2*(-x4)*x3*y3*eta;
             1.0/2.0*z1*zet*x1*(-y3)*y1             -1.0/2.0*z1*zet*x3*xi*eta   1.0/2.0*z2*x3*xi*(-y3)*y1;
                       -1.0/2.0*z3*zet*xi*y1*eta   1.0/2.0*z3*zet*(-x4)*x3*y2   1.0/2.0*z4*(-x4)*x3*y1*eta;
             1.0/2.0*z3*zet*x2*(-y3)*y1             -1.0/2.0*z3*zet*x4*xi*eta   1.0/2.0*z4*x4*xi*(-y3)*y1;
                       -1.0/2.0*z3*zet*xi*y3*eta   1.0/2.0*z3*zet*(-x4)*x3*y4   1.0/2.0*z4*(-x4)*x3*y3*eta;
             1.0/2.0*z3*zet*x1*(-y3)*y1             -1.0/2.0*z3*zet*x3*xi*eta   1.0/2.0*z4*x3*xi*(-y3)*y1;
             1.0/2.0*(-z3)*z1*x1*y1*eta            1.0/2.0*(-z3)*z1*x3*xi*y2             -1.0/2.0*zet*x3*xi*y1*eta;
             1.0/2.0*(-z3)*z1*x2*y1*eta            1.0/2.0*(-z3)*z1*x4*xi*y2             -1.0/2.0*zet*x4*xi*y1*eta;
             1.0/2.0*(-z3)*z1*x2*y3*eta            1.0/2.0*(-z3)*z1*x4*xi*y4             -1.0/2.0*zet*x4*xi*y3*eta;
             1.0/2.0*(-z3)*z1*x1*y3*eta            1.0/2.0*(-z3)*z1*x3*xi*y4             -1.0/2.0*zet*x3*xi*y3*eta;
                      -z1*zet*xi*(-y3)*y1            -z1*zet*(-x4)*x3*eta  z2*(-x4)*x3*(-y3)*y1;
                      -(-z3)*z1*xi*y1*eta  (-z3)*z1*(-x4)*x3*y2            -zet*(-x4)*x3*y1*eta;
            (-z3)*z1*x2*(-y3)*y1            -(-z3)*z1*x4*xi*eta            -zet*x4*xi*(-y3)*y1;
                      -(-z3)*z1*xi*y3*eta  (-z3)*z1*(-x4)*x3*y4            -zet*(-x4)*x3*y3*eta;
            (-z3)*z1*x1*(-y3)*y1            -(-z3)*z1*x3*xi*eta            -zet*x3*xi*(-y3)*y1;
                      -z3*zet*xi*(-y3)*y1            -z3*zet*(-x4)*x3*eta  z4*(-x4)*x3*(-y3)*y1;
               -2.0*(-z3)*z1*xi*(-y3)*y1     -2.0*(-z3)*z1*(-x4)*x3*eta     -2.0*zet*(-x4)*x3*(-y3)*y1;
           ];
    return reshape(val, 27, 3);
end

function privboundaryconn(self::FESetH27)
    return  [self.conn[:, [1,  4,  3,  2,  12,  11,  10,  9,  21]];
             self.conn[:, [1,  2,  6,  5,  9,  18,  13,  17,  22]];
             self.conn[:, [2,  3,  7,  6,  10,  19,  14,  18,  23]];
             self.conn[:, [3,  4,  8,  7,   11,  20,  15,  19,  24]];
             self.conn[:, [4,  1,  5,  8,  12,  17,  16,  20,  25]];
             self.conn[:, [6,  7,  8,  5,  14,  15,  16,  13,  26]]];
end

function privboundaryfe(self::FESetH27)
    return FESetQ9;
end

"""
    FESetT4

Type for sets of volume-like of tetrahedral finite elements with four nodes.
"""
mutable struct FESetT4 <: FESet3Manifold
  @add_FESet_fields

  function    FESetT4(conn::FIntMat=[])
    nodesperelem::FInt  = 4
    @assert (size(conn, 2) == nodesperelem) "Number of nodes per element mismatched"
    # Need to make a COPY of the input arrays
    self =new(nodesperelem, deepcopy(conn), deepcopy(FInt[]))
    setlabel!(self, 0)
    return self
  end
end

function privbfun(self::FESetT4,  param_coords::FFltVec)
    # Evaluate the basis function matrix for an 3-node triangle.
     val = [(1 - param_coords[1] - param_coords[2] - param_coords[3]);
            param_coords[1];
            param_coords[2];
            param_coords[3]];
    return reshape(val, 4, 1)
end

function privbfundpar(self::FESetT4,  param_coords::FFltVec)
    # Evaluate the derivatives of the basis function matrix.
    val = [-1. -1. -1.;
           +1.  0.  0.;
           0. +1.  0.;
           0.  0. +1.];
    return reshape(val, 4, 3)
end

function privboundaryconn(self::FESetT4)
    return [self.conn[:, [1,  3,  2]];self.conn[:, [1,  2,  4]];self.conn[:, [2,  3,  4]];self.conn[:, [1,  4,  3]]];
end

function privboundaryfe(self::FESetT4)
    return FESetT3;
end

"""
    FESetT10

Type for sets of volume-like of tetrahedral finite elements with 10 nodes.
"""
mutable struct FESetT10 <: FESet3Manifold
  @add_FESet_fields

  function FESetT10(conn::FIntMat=[])
    nodesperelem::FInt  = 10
    @assert (size(conn, 2) == nodesperelem) "Number of nodes per element mismatched"
    # Need to make a COPY of the input arrays
    self =new(nodesperelem, deepcopy(conn), deepcopy(FInt[]))
    setlabel!(self, 0)
    return self
  end
end

function privbfun(self::FESetT10,  param_coords::FFltVec)
    # Evaluate the basis function matrix for an 3-node triangle.
    r = param_coords[1];
    s = param_coords[2];
    t = param_coords[3];
    val = [(1-r-s-t) * (2*(1-r-s-t)-1);
           r * (2*r-1);
           s * (2*s-1);
           t * (2*t-1);
           4*(1-r-s-t)*r;
           4*r*s;
           4*s*(1-r-s-t);
           4*(1-r-s-t)*t;
           4*r*t;
           4*s*t;
           ];
    return reshape(val, 10, 1)
end

function privbfundpar(self::FESetT10,  param_coords::FFltVec)
    # Evaluate the derivatives of the basis function matrix.
    r = param_coords[1];
    s = param_coords[2];
    t = param_coords[3];
            val = [-3+4*r+4*s+4*t   4*r-1 0   0  -8*r+4-4*s-4*t  4*s -4*s  -4*t  4*t  0;
                -3+4*r+4*s+4*t  0  4*s-1  0   -4*r   4*r 4-4*r-8*s-4*t -4*t  0 4*t;
                -3+4*r+4*s+4*t  0  0  4*t-1   -4*r   0  -4*s -8*t+4-4*r-4*s  4*r  4*s
                ]';
    return reshape(val, 10, 3)
end

function privboundaryconn(self::FESetT10)
    return [self.conn[:, [1,  3,  2,  7,  6,  5]];
    self.conn[:, [1,  2,  4,  5,  9,  8]];
    self.conn[:, [2,  3,  4,  6,  10,  9]];
    self.conn[:, [3,  1,  4,  7,  8,  10]]];
end

function privboundaryfe(self::FESetT10)
    return FESetT6;
end

end
