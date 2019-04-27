"""
    FESetModule

Module for finite element sets.
"""
module FESetModule

import Base.count
import Base.cat
import LinearAlgebra: mul!, Transpose
At_mul_B!(C, A, B) = mul!(C, Transpose(A), B)
import LinearAlgebra: norm, cross
using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict

"""
    AbstractFESet{NODESPERELEM}

Abstract type of a finite element set. Parameterized with the number of of the nodes per element.
"""
abstract type AbstractFESet{NODESPERELEM} end
"""
    AbstractFESet0Manifold{NODESPERELEM} <: FESet{NODESPERELEM}

Abstract type of a finite element set for 0-dimensional manifolds (points).
Parameterized with the number of of the nodes per element.
"""
abstract type AbstractFESet0Manifold{NODESPERELEM} <: AbstractFESet{NODESPERELEM} end
"""
    AbstractFESet1Manifold{NODESPERELEM} <: FESet{NODESPERELEM}

Abstract type of a finite element set for 1-dimensional manifolds (curves).
Parameterized with the number of of the nodes per element.
"""
abstract type AbstractFESet1Manifold{NODESPERELEM} <: AbstractFESet{NODESPERELEM} end
"""
    AbstractFESet2Manifold{NODESPERELEM} <: FESet{NODESPERELEM}

Abstract type of a finite element set for 2-dimensional manifolds (surfaces).
Parameterized with the number of of the nodes per element.
"""
abstract type AbstractFESet2Manifold{NODESPERELEM} <: AbstractFESet{NODESPERELEM} end
"""
    AbstractFESet3Manifold{NODESPERELEM} <: FESet{NODESPERELEM}

Abstract type of a finite element set for 3-dimensional manifolds (solids).
Parameterized with the number of of the nodes per element.
"""
abstract type AbstractFESet3Manifold{NODESPERELEM} <: AbstractFESet{NODESPERELEM} end

"""
    add_FESet_fields(NODESPERELEM)

Generate standard fields for the finite element set.
"""
macro add_FESet_fields(NODESPERELEM)
    return esc(:(
    conn::Array{NTuple{$NODESPERELEM, FInt}, 1};
    label::FIntVec; )
    )
end
# show(macroexpand(:(@add_FESet_fields 6)))

"""
    define_FESet(NAME, MANIFOLD, NODESPERELEM)

Define the concrete type for a finite element set.
"""
macro define_FESet(NAME, MANIFOLD, NODESPERELEM)
    return esc(:(
        mutable struct $NAME <: $MANIFOLD{$NODESPERELEM}
            @add_FESet_fields $NODESPERELEM
            function $NAME(conn::FIntMat)
                self = new(NTuple{$NODESPERELEM, FInt}[], FInt[])
                self = fromarray!(self, conn)
                setlabel!(self, 0)
                return self
            end
        end
    ))
end
# show(macroexpand(:(@define_FESet FESetT10 AbstractFESet3Manifold 10)))

"""
    nodesperelem(fes::AbstractFESet{NODESPERELEM}) where {NODESPERELEM}

Provide the number of nodes per element.  
"""
nodesperelem(fes::AbstractFESet{NODESPERELEM}) where {NODESPERELEM} = NODESPERELEM

"""
    manifdim(me)

Get the manifold dimension.
"""
manifdim(me::AbstractFESet0Manifold{NODESPERELEM}) where {NODESPERELEM} = 0
manifdim(me::AbstractFESet1Manifold{NODESPERELEM}) where {NODESPERELEM} = 1
manifdim(me::AbstractFESet2Manifold{NODESPERELEM}) where {NODESPERELEM} = 2
manifdim(me::AbstractFESet3Manifold{NODESPERELEM}) where {NODESPERELEM} = 3

"""
    count(self::T)::FInt where {T<:AbstractFESet}

Get the number of individual connectivities in the FE set.
"""
count(self::T) where {T<:AbstractFESet} = length(self.conn)

"""
    fromarray!(self::AbstractFESet{NODESPERELEM}, conn::FIntMat) where {NODESPERELEM}

Set  the connectivity from an integer array.  
"""
function fromarray!(self::AbstractFESet{NODESPERELEM}, conn::FIntMat) where {NODESPERELEM}
    @assert size(conn, 2) == NODESPERELEM
    self.conn = fill(tuple(fill(0, NODESPERELEM)...), size(conn, 1));
    for i = 1:length(self.conn)
        self.conn[i] = ntuple(y -> conn[i, y], NODESPERELEM);
    end
    return self
end

"""
    connasarray(self::AbstractFESet{NODESPERELEM}) where {NODESPERELEM}

Return the connectivity as an integer array (matrix), where the number of rows
matches the number of connectivities in the set.
"""
function connasarray(self::AbstractFESet{NODESPERELEM}) where {NODESPERELEM}
    conn = fill(zero(FInt), length(self.conn), NODESPERELEM)
    for i = 1:size(conn, 1)
        for j = 1:NODESPERELEM
            conn[i, j] = self.conn[i][j]
        end
    end
    return conn
end

"""
    boundaryconn(self::T) where {T<:AbstractFESet}

Get boundary connectivity.
"""
function boundaryconn(self::T) where {T<:AbstractFESet}
    return privboundaryconn(self)
end

"""
    boundaryfe(self::T) where {T<:AbstractFESet}

Return the constructor of the type of the boundary finite element.
"""
function boundaryfe(self::T) where {T<:AbstractFESet}
    return privboundaryfe(self)
end

"""
    bfun(self::T,  param_coords::FFltVec)::FFltMat where {T<:AbstractFESet}

Compute the values of the basis functions at a given parametric coordinate.
One basis function per row.
"""
function bfun(self::T,  param_coords::FFltVec)::FFltMat where {T<:AbstractFESet}
    return privbfun(self, param_coords)
end

"""
    bfundpar(self::T,  param_coords::FFltVec)::FFltMat where {T<:AbstractFESet}

Compute the values of the basis function gradients with respect to the
parametric coordinates at a given parametric coordinate. One basis function
gradients per row.
"""
function bfundpar(self::T,  param_coords::FFltVec)::FFltMat where {T<:AbstractFESet}
    return privbfundpar(self, param_coords)
end

"""
    setlabel!(self::T, val::FInt) where {T<:AbstractFESet}

Set the label of the entire finite elements set.
"""
function setlabel!(self::T, val::FInt) where {T<:AbstractFESet}
    self.label = zeros(FInt, size(self.conn, 1));
    fill!(self.label, val)
    return self
end

"""
    setlabel!(self::T, val::FIntVec) where {T<:AbstractFESet}

Set the labels of individual elements.
"""
function setlabel!(self::T, val::FIntVec) where {T<:AbstractFESet}
    #    Set the label of this set.
    @assert size(self.conn, 1)==length(val) "Must get one  label per finite element connectivity"
    self.label=zeros(FInt, size(self.conn, 1));
    copyto!(self.label, val);
    return self
end

"""
    subset(self::T, L::FIntVec) where {T<:AbstractFESet}

Extract a subset of the finite elements from the given finite element set.
"""
function subset(self::T, L::FIntVec) where {T<:AbstractFESet}
    result = deepcopy(self)
    result.conn = deepcopy(self.conn[L])
    result.label = deepcopy(self.label[L])
    return  result
end

"""
    cat(self::T,  other::T) where {T<:AbstractFESet}

Concatenate the connectivities of two FE sets.
"""
function cat(self::T,  other::T) where {T<:AbstractFESet}
    @assert nodesperelem(self) == nodesperelem(other)
    result =deepcopy(self)
    result.conn = vcat(self.conn, other.conn);
    setlabel!(result, vcat(self.label, other.label))
    return result
end

"""
    updateconn!(self::T, newids::FIntVec) where {T<:AbstractFESet}

Update the connectivity after the IDs of nodes changed.

`newids`= new node IDs. Note that indexes in the conn array "point"
_into_ the  `newids` array. After the connectivity was updated
this will no longer be true!
"""
function updateconn!(self::T, newids::FIntVec) where {T<:AbstractFESet}
    conn = connasarray(self)
    for i=1:size(conn, 1)
        for j=1:size(conn, 2)
            conn[i, j]=newids[conn[i, j]];
        end
    end
    return fromarray!(self, conn)
end

"""
    inparametric(self::AbstractFESet, param_coords::FFltVec)

Are given parametric coordinates inside the element parametric domain?

Return a Boolean: is the point inside, true or false?
"""
function inparametric(self::T, param_coords::FFltVec; tolerance = 0.0) where {T<:AbstractFESet}
    return privinparametric(self, param_coords, tolerance)
end

"""
    map2parametric(self::T, x::FFltMat, pt::FFltVec;
        tolerance = 0.001, maxiter =5) where {T<:AbstractFESet}

Map a spatial location to parametric coordinates.

- `x`=array of spatial coordinates of the nodes, size(x) = nbfuns x dim,
- `c`= spatial location
- `tolerance` = tolerance in parametric coordinates; default is 0.001.

# Return
- `success` = Boolean flag, true if successful, false otherwise.
- `pc` = Returns a row array of parametric coordinates if the solution was
  successful, otherwise NaN are returned.
"""
function map2parametric(self::T, x::FFltMat, pt::FFltVec; tolerance = 0.001, maxiter =5) where {T<:AbstractFESet}
    sdim = size(x, 2); # number of space dimensions
    mdim = manifdim(self); # manifold dimension of the element
    ppc = zeros(mdim);
    pc = deepcopy(ppc);
    J = [i==j ? one(FFlt) : zero(FFlt) for i=1:sdim, j=1:mdim] # Jac. matx: "identity"
    loc = fill(zero(FFlt), (1, sdim)) # point location buffer
    success = true
    for i = 1:maxiter
        gradNparams = privbfundpar(self, pc);
        At_mul_B!(J, x, gradNparams); # calculate the Jacobian matrix
        N = privbfun(self, pc);
        At_mul_B!(loc, N, x);# Iterated point location
        pc[:] = ppc .- (J\(vec(loc) - pt))
        if (norm(pc-ppc) < tolerance)
            return pc, success;
        end
        ppc[:] = pc[:];
    end
    success = false
    return pc, success;
end

"""
    centroidparametric(self::T) where {T<:AbstractFESet}

Return the parametric coordinates  of the centroid of the element.
"""
function centroidparametric(self::T) where {T<:AbstractFESet}
    return privcentroidparametric(self)
end

"""
    Jacobian(self::T, J::FFltMat)::FFlt where {T<:AbstractFESet1Manifold}

Evaluate the curve Jacobian.

- `J` = Jacobian matrix, columns are tangent to parametric coordinates curves.
"""
function Jacobian(self::T, J::FFltMat)::FFlt where {T<:AbstractFESet0Manifold}
    return 1.0::FFlt;
end

"""
    Jacobian(self::T, J::FFltMat)::FFlt where {T<:AbstractFESet1Manifold}

Evaluate the curve Jacobian.

- `J` = Jacobian matrix, columns are tangent to parametric coordinates curves.
"""
function Jacobian(self::T, J::FFltMat)::FFlt where {T<:AbstractFESet1Manifold}
    sdim,  ntan = size(J);
    @assert ntan == 1 "Expected number of tangent vectors: 1"
    return norm(J)::FFlt;
end

"""
    gradN!(self::AbstractFESet1Manifold, gradN::FFltMat, gradNparams::FFltMat,
      redJ::FFltMat)

Compute the gradient of the basis functions with the respect to
the "reduced" spatial coordinates.

- `gradN`= output,  matrix of gradients,  one per row
- `gradNparams`= matrix of gradients with respect to parametric coordinates, one per row
- `redJ`= reduced Jacobian matrix `redJ=transpose(Rm)*J`
"""
function gradN!(self::AbstractFESet1Manifold, gradN::FFltMat, gradNparams::FFltMat, redJ::FFltMat)
    for r=1:size(gradN, 1)
        gradN[r, 1]=  gradNparams[r, 1]/redJ[1, 1]
    end
end

"""
    Jacobian(self::T, J::FFltMat)::FFlt where {T<:AbstractFESet2Manifold}

Evaluate the curve Jacobian.

- `J` = Jacobian matrix, columns are tangent to parametric coordinates curves.
"""
function Jacobian(self::T, J::FFltMat)::FFlt where {T<:AbstractFESet2Manifold}
    sdim,  ntan = size(J);
    @assert ntan == 2 "Expected number of tangent vectors: 2"
    if sdim == ntan
        @inbounds Jac = (J[1, 1]*J[2, 2] - J[2, 1]*J[1, 2])
        return Jac::FFlt;# is det(J);% Compute the Jacobian
    else
        return norm(cross(J[:, 1], J[:, 2]))::FFlt;
    end
end

"""
    gradN!(self::AbstractFESet2Manifold, gradN::FFltMat, gradNparams::FFltMat,
      redJ::FFltMat)

Compute the gradient of the basis functions with the respect to
the "reduced" spatial coordinates.

- `gradN`= output,  matrix of gradients,  one per row
- `gradNparams`= matrix of gradients with respect to parametric coordinates,
  one per row
- `redJ`= reduced Jacobian matrix `redJ=transpose(Rm)*J`
"""
function gradN!(self::AbstractFESet2Manifold, gradN::FFltMat, gradNparams::FFltMat, redJ::FFltMat)
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
    Jacobian(self::T, J::FFltMat)::FFlt where {T<:AbstractFESet3Manifold}

Evaluate the volume Jacobian.

`J` = Jacobian matrix, columns are tangent to parametric coordinates curves.
"""
function Jacobian(self::T, J::FFltMat)::FFlt where {T<:AbstractFESet3Manifold}
    sdim,  ntan = size(J);
    @assert (ntan == 3) && (sdim == 3) "Expected number of tangent vectors: 3"
    #Jac = det(J);# Compute the Jacobian
    # The unrolled version
    return (+J[1, 1]*(J[2, 2]*J[3, 3]-J[3, 2]*J[2, 3])
    -J[1, 2]*(J[2, 1]*J[3, 3]-J[2, 3]*J[3, 1])
    +J[1, 3]*(J[2, 1]*J[3, 2]-J[2, 2]*J[3, 1]) )::FFlt;
end

"""
    gradN!(self::AbstractFESet3Manifold, gradN::FFltMat, gradNparams::FFltMat,
      redJ::FFltMat)

Compute the gradient of the basis functions with the respect to
the "reduced" spatial coordinates.

- `gradN`= output,  matrix of gradients,  one per row
- `gradNparams`= matrix of gradients with respect to parametric coordinates,
  one per row
- `redJ`= reduced Jacobian matrix `redJ=transpose(Rm)*J`
"""
function gradN!(self::AbstractFESet3Manifold, gradN::FFltMat, gradNparams::FFltMat, redJ::FFltMat)
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


    
################################################################################
################################################################################

"""
    FESetP1

Type for sets of point-like of finite elements.
"""
@define_FESet FESetP1 AbstractFESet0Manifold 1

privbfun(self::FESetP1,  param_coords::FFltVec) = reshape([1.0], 1, 1) # make sure this is a matrix
privbfundpar(self::FESetP1,  param_coords::FFltVec) = zeros(1, 0)


function privinparametric(self::FESetP1, param_coords::FFltVec, tolerance::FFlt)
    return  param_coords[1] == 0.0;
end

function privcentroidparametric(self::FESetP1)
    return vec([0.0])
end

function map2parametric(self::FESetP1, x::FFltMat, pt::FFltVec;
    tolerance = 0.001, maxiter =5)
    success = false; pc = [0.0]
    if norm(vec(x) - pt) < tolerance
        success = true
    end
    return pc, success;
end

################################################################################
################################################################################

"""
    FESetL2

Type for sets of curve-like finite elements with two nodes.
"""
@define_FESet FESetL2 AbstractFESet1Manifold 2

privbfun(self::FESetL2,  param_coords::FFltVec) = reshape([(1. - param_coords[1]); (1. + param_coords[1])] / 2.0, 2, 1) # make sure this is a matrix
privbfundpar(self::FESetL2,  param_coords::FFltVec) = reshape([-1.0; +1.0]/2.0, 2, 1)

function privboundaryconn(self::FESetL2)
    conn = connasarray(self)
    return [conn[:, 1]; conn[:, 2]];
end

function privboundaryfe(self::FESetL2)
    return FESetP1;
end

function privinparametric(self::FESetL2, param_coords::FFltVec, tolerance::FFlt)
    s = findall(v -> (v >= (-1.0 - tolerance)) && (v <= (+1.0 + tolerance)), param_coords);
    return  (length(s) == length(param_coords));
end

function privcentroidparametric(self::FESetL2)
    return vec([0.0])
end

################################################################################
################################################################################

"""
    FESetL3

Type for sets of curve-like of finite elements with three nodes.
"""
@define_FESet FESetL3 AbstractFESet1Manifold 3

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
    conn = connasarray(self)
    return [conn[:, 1]; conn[:, 2]];
end

function privboundaryfe(self::FESetL3)
    return FESetP1;
end

function privinparametric(self::FESetL3, param_coords::FFltVec, tolerance::FFlt)
    s = findall(v -> (v >= (-1.0 - tolerance)) && (v <= (+1.0 + tolerance)), param_coords);
    return  (length(s) == length(param_coords));
end

function privcentroidparametric(self::FESetL3)
    return vec([0.0])
end

################################################################################
################################################################################

"""
    FESetT3

Type for sets of surface-like triangular finite elements with three nodes.
"""
@define_FESet FESetT3 AbstractFESet2Manifold 3

function privbfun(self::FESetT3,  param_coords::FFltVec)
    # Evaluate the basis function matrix for an 3-node triangle.
    return reshape([(1 - param_coords[1] - param_coords[2]); param_coords[1]; param_coords[2]], 3, 1); # Make sure this is a matrix
end

function privbfundpar(self::FESetT3,  param_coords::FFltVec)
    # Evaluate the derivatives of the basis function matrix.
    return  [-1. -1.;  +1.  0.;  0. +1.];
end

function privboundaryconn(self::FESetT3)
    conn = connasarray(self)
    return [conn[:, 1:2]; conn[:, 2:3]; conn[:, [3, 1]]];
end

function privboundaryfe(self::FESetT3)
    # Get  the constructor of the class of the  boundary finite element.
    return FESetL2;
end

function privinparametric(self::FESetT3, param_coords::FFltVec, tolerance::FFlt)
    s = findall(v -> (v >= (-tolerance)) && (v <= (+1.0 + tolerance)), param_coords);
    return  (length(s) == length(param_coords)) &&
        (-3.0 * tolerance <= sum(param_coords) <= 1.0 + 3.0 * tolerance);
end

function privcentroidparametric(self::FESetT3)
    return vec([1.0/3 1.0/3])
end

################################################################################
################################################################################

"""
    FESetQ4

Type for sets of surface-like quadrilateral finite elements with four nodes.
"""
@define_FESet FESetQ4 AbstractFESet2Manifold 4

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
    conn = connasarray(self)
    return [conn[:, 1:2]; conn[:, 2:3]; conn[:, 3:4]; conn[:, [4, 1]]];
end

function privboundaryfe(self::FESetQ4)
    # Get  the constructor of the class of the  boundary finite element.
    return FESetL2;
end

function privinparametric(self::FESetQ4, param_coords::FFltVec, tolerance::FFlt)
    s = findall(v -> (v >= (-1.0 - tolerance)) && (v <= (+1.0 + tolerance)), param_coords);
    return  (length(s) == length(param_coords));
end

function privcentroidparametric(self::FESetQ4)
    return vec([0.0 0.0])
end

################################################################################
################################################################################

"""
    FESetQ9

Type for sets of surface-like quadrilateral finite elements with nine nodes.
"""
@define_FESet FESetQ9 AbstractFESet2Manifold 9

function privbfun(self::FESetQ9,  param_coords::FFltVec)
    # Evaluate the basis function matrix for an 9-node quadrilateral.
    xi=param_coords[1];
    xis = [(xi-1)*xi/2;  (xi+1)*xi/2;  -(xi+1)*(xi-1)];
    eta=param_coords[2];
    etas = [(eta-1)*eta/2;  (eta+1)*eta/2;  -(eta+1)*(eta-1)];
    xisetas=(xis*etas');
    val = xisetas[[1     2     5     4     3     8     6     7     9]];
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
    val =  [dxisetas[ix]'             xisdetas[ix]'];
    return reshape(val, 9, 2); # Make sure this is a matrix
end

function privboundaryconn(self::FESetQ9)
    conn = connasarray(self)
    return [conn[:, [1, 2, 5]]; conn[:, [2, 3, 6]]; conn[:, [3, 4, 7]]; conn[:, [4, 1, 8]]];
end

function privboundaryfe(self::FESetQ9)
    return FESetL3;
end

function privinparametric(self::FESetQ9, param_coords::FFltVec, tolerance::FFlt)
    s = findall(v -> (v >= (-1.0 - tolerance)) && (v <= (+1.0 + tolerance)), param_coords);
    return  (length(s) == length(param_coords));
end

function privcentroidparametric(self::FESetQ9)
    return vec([0.0 0.0])
end

################################################################################
################################################################################

"""
    FESetQ8

Type for sets of surface-like quadrilateral finite elements with eight nodes.
"""
@define_FESet FESetQ8 AbstractFESet2Manifold 8

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
    conn = connasarray(self)
    return [conn[:, [1, 2, 5]]; conn[:, [2, 3, 6]]; conn[:, [3, 4, 7]]; conn[:, [4, 1, 8]]];
end

function privboundaryfe(self::FESetQ8)
    return FESetL3;
end

function privinparametric(self::FESetQ8, param_coords::FFltVec, tolerance::FFlt)
    s = findall(v -> (v >= (-1.0 - tolerance)) && (v <= (+1.0 + tolerance)), param_coords);
    return  (length(s) == length(param_coords));
end

function privcentroidparametric(self::FESetQ8)
    return vec([0.0 0.0])
end

################################################################################
################################################################################

"""
    FESetT6

Type for sets of surface-like triangular finite elements with six nodes.
"""
@define_FESet FESetT6 AbstractFESet2Manifold 6

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
    conn = connasarray(self)
    return [conn[:, [1,  2,  4]]; conn[:, [2,  3,  5]]; conn[:, [3,  1,  6]]];
end

function privboundaryfe(self::FESetT6)
    # Get  the constructor of the class of the  boundary finite element.
    return FESetL3;
end

function privinparametric(self::FESetT6, param_coords::FFltVec, tolerance::FFlt)
    s = findall(v -> (v >= (-tolerance)) && (v <= (+1.0 + tolerance)), param_coords);
    return  (length(s) == length(param_coords)) &&
      (-3.0 * tolerance <= sum(param_coords) <= 1.0 + 3.0 * tolerance);
end

function privcentroidparametric(self::FESetT6)
    return vec([1.0/3 1.0/3])
end

################################################################################
################################################################################

"""
    FESetH8

Type for sets of volume-like hexahedral finite elements with eight nodes.
"""
@define_FESet FESetH8 AbstractFESet3Manifold 8

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
    conn = connasarray(self)
    return  [conn[:, [1,  4,  3,  2]];
             conn[:, [1,  2,  6,  5]];
             conn[:, [2,  3,  7,  6]];
             conn[:, [3,  4,  8,  7]];
             conn[:, [4,  1,  5,  8]];
             conn[:, [6,  7,  8,  5]]];
end

function privboundaryfe(self::FESetH8)
    # Get  the constructor of the class of the  boundary finite element.
    return FESetQ4;
end

function privinparametric(self::FESetH8, param_coords::FFltVec, tolerance::FFlt)
    s = findall(v -> (v >= (-1.0 - tolerance)) && (v <= (+1.0 + tolerance)), param_coords);
    return  (length(s) == length(param_coords));
end

function privcentroidparametric(self::FESetH8)
    return vec([0.0 0.0 0.0])
end

################################################################################
################################################################################

"""
    FESetH20

Type for sets of volume-like hexahedral finite elements with 20 nodes.
"""
@define_FESet FESetH20 AbstractFESet3Manifold 20

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
    conn = connasarray(self)
    return  [conn[:, [1,  4,  3,  2,  12,  11,  10,  9]];
             conn[:, [1,  2,  6,  5,  9,  18,  13,  17]];
             conn[:, [2,  3,  7,  6,  10,  19,  14,  18]];
             conn[:, [3,  4,  8,  7,   11,  20,  15,  19]];
             conn[:, [4,  1,  5,  8,  12,  17,  16,  20]];
             conn[:, [6,  7,  8,  5,  14,  15,  16,  13]]];
end

function privboundaryfe(self::FESetH20)
    # Get  the constructor of the class of the  boundary finite element.
    return FESetQ8;
end

function privinparametric(self::FESetH20, param_coords::FFltVec, tolerance::FFlt)
    s = findall(v -> (v >= (-1.0 - tolerance)) && (v <= (+1.0 + tolerance)), param_coords);
    return  (length(s) == length(param_coords));
end

function privcentroidparametric(self::FESetH20)
    return vec([0.0 0.0 0.0])
end

################################################################################
################################################################################

"""
    FESetH27

Type for sets of volume-like hexahedral finite elements with 27 nodes.
"""
@define_FESet FESetH27 AbstractFESet3Manifold 27

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
    conn = connasarray(self)
    return  [conn[:, [1,  4,  3,  2,  12,  11,  10,  9,  21]];
             conn[:, [1,  2,  6,  5,  9,  18,  13,  17,  22]];
             conn[:, [2,  3,  7,  6,  10,  19,  14,  18,  23]];
             conn[:, [3,  4,  8,  7,   11,  20,  15,  19,  24]];
             conn[:, [4,  1,  5,  8,  12,  17,  16,  20,  25]];
             conn[:, [6,  7,  8,  5,  14,  15,  16,  13,  26]]];
end

function privboundaryfe(self::FESetH27)
    return FESetQ9;
end

function privinparametric(self::FESetH27, param_coords::FFltVec, tolerance::FFlt)
    s = findall(v -> (v >= (-1.0 - tolerance)) && (v <= (+1.0 + tolerance)), param_coords);
    return  (length(s) == length(param_coords));
end

function privcentroidparametric(self::FESetH27)
    return vec([0.0 0.0 0.0])
end

################################################################################
################################################################################

"""
    FESetT4

Type for sets of volume-like tetrahedral finite elements with four nodes.
"""
@define_FESet FESetT4 AbstractFESet3Manifold 4

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
    conn = connasarray(self)
    return [conn[:, [1,  3,  2]]; conn[:, [1,  2,  4]]; conn[:, [2,  3,  4]]; conn[:, [1,  4,  3]]];
end

function privboundaryfe(self::FESetT4)
    return FESetT3;
end

function privinparametric(self::FESetT4, param_coords::FFltVec, tolerance::FFlt)
    s = findall(v -> (v >= (-tolerance)) && (v <= (+1.0 + tolerance)), param_coords);
    return  (length(s) == length(param_coords)) &&
      (-3.0 * tolerance <= sum(param_coords) <= 1.0 + 3.0 * tolerance);
end

function privcentroidparametric(self::FESetT4)
    return vec([1.0/4 1.0/4 1.0/4])
end

################################################################################
################################################################################

"""
    FESetT10

Type for sets of volume-like tetrahedral finite elements with 10 nodes.
"""
@define_FESet FESetT10 AbstractFESet3Manifold 10

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
    conn = connasarray(self)
    return [conn[:, [1,  3,  2,  7,  6,  5]];
            conn[:, [1,  2,  4,  5,  9,  8]];
            conn[:, [2,  3,  4,  6,  10,  9]];
            conn[:, [3,  1,  4,  7,  8,  10]]];
end

function privboundaryfe(self::FESetT10)
    return FESetT6;
end

function privinparametric(self::FESetT10, param_coords::FFltVec, tolerance::FFlt)
    s = findall(v -> (v >= (-tolerance)) && (v <= (+1.0 + tolerance)), param_coords);
    return  (length(s) == length(param_coords)) &&
      (-3.0 * tolerance <= sum(param_coords) <= 1.0 + 3.0 * tolerance);
end

function privcentroidparametric(self::FESetT10)
    return vec([1.0/4 1.0/4 1.0/4])
end

################################################################################
################################################################################

end
