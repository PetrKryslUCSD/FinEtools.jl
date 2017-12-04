"""
    IntegDataModule

Module to manage geometry data.
"""
module IntegDataModule

using FinEtools.FTypesModule
import FinEtools.FESetModule: FESet, FESet0Manifold, FESet1Manifold, FESet2Manifold, FESet3Manifold, Jacobian, bfun, bfundpar
import FinEtools.FENodeSetModule: FENodeSet
import FinEtools.IntegRuleModule: IntegRule

"""
    IntegData{T<:FESet}

Integration data.

T = type of finite element set.  The type of the FE set will be dependent upon the
operations required. For instance, for interior (volume) integrals such as body
load or the stiffness hexahedral H8 may be used, whereas for boundary  (surface)
integrals quadrilateral Q4 would be needed.
"""
mutable struct IntegData{S<:FESet, F<:Function}
    fes::S # finite element set object
    integration_rule::IntegRule  # integration rule object
    otherdimension::F
    axisymmetric::Bool
end

"""
    IntegData(fes::S, integration_rule::IntegRule) where {S<:FESet}

Construct with the default orientation matrix (identity), and the other
dimension  being the default 1.0.
"""
function  IntegData(fes::S, integration_rule::IntegRule) where {S<:FESet}
    return IntegData(fes, integration_rule, otherdimensionunity, false)
end

"""
    IntegData(fes::S, integration_rule::IntegRule,
      otherdimension::FFlt) where {S<:FESet}

Construct with the default orientation matrix (identity), and constant other
dimension.
"""
function  IntegData(fes::S, integration_rule::IntegRule, otherdimension::FFlt) where {S<:FESet}
    function otherdimensionfu(loc::FFltMat,
        conn::CC, N::FFltMat)::FFlt where {CC}
        return otherdimension::FFlt
    end
    return IntegData(fes, integration_rule, otherdimensionfu, false)
end

"""
    IntegData(fes::S, integration_rule::IntegRule,
      axisymmetric::Bool) where {S<:FESet}

Construct with the default orientation matrix (identity), for axially symmetric
models. The other dimension is  the default  unity (1.0).
"""
function IntegData(fes::S, integration_rule::IntegRule, axisymmetric::Bool) where {S<:FESet}
    return IntegData(fes, integration_rule, otherdimensionunity, axisymmetric)
end

"""
    IntegData(fes::S, integration_rule::IntegRule, axisymmetric::Bool,
      otherdimension::FFlt) where {S<:FESet}

Construct for axially symmetric models. The other dimension is given as a number.
"""
function IntegData(fes::S, integration_rule::IntegRule, axisymmetric::Bool, otherdimension::FFlt) where {S<:FESet}
    function otherdimensionfu(loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {CC}
        return otherdimension::FFlt
    end
    return IntegData(fes, integration_rule, otherdimensionfu, axisymmetric)
end

"""
    otherdimensionunity(loc::FFltMat, conn::CC, N::FFltMat)::FFlt
              where {CC}

Evaluate the other dimension: default is 1.0.
"""
function otherdimensionunity(loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {CC}
    return 1.0
end

"""
    Jacobianpoint(self::IntegData{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                N::FFltMat)::FFlt where {T<:FESet0Manifold, CC}

Evaluate the point Jacobian.

`J` = Jacobian matrix
`loc` = location of the quadrature point in physical coordinates,
`conn` = connectivity of the element,
`N` = matrix of basis function values at the quadrature point.
"""
function Jacobianpoint(self::IntegData{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {T<:FESet0Manifold, CC}
    return Jacobian(self.fes, J)::FFlt
end

"""
    Jacobiancurve(self::IntegData{T}, J::FFltMat,
        loc::FFltMat, conn::CC,
        N::FFltMat)::FFlt where {T<:FESet0Manifold, CC}

Evaluate the curve Jacobian.

`J` = Jacobian matrix
`loc` = location of the quadrature point in physical coordinates,
`conn` = connectivity of the element,
`N` = matrix of basis function values at the quadrature point.
"""
function Jacobiancurve(self::IntegData{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {T<:FESet0Manifold, CC}
    Jac = Jacobianpoint(self, J, loc, conn, N)
    if self.axisymmetric
        return Jac*2*pi*loc[1];
    else
        return Jac*self.otherdimension(loc, conn,  N)
    end
end

"""
    Jacobiansurface(self::IntegData{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                N::FFltMat)::FFlt where {T<:FESet0Manifold, CC}

Evaluate the surface Jacobian.

For the zero-dimensional cell,  the surface Jacobian is
    (i) the product of the point Jacobian and the other dimension
    (units of length squared);
    or,  when used as axially symmetric
    (ii) the product of the point Jacobian and the circumference of
    the circle through the point `loc` times the other dimension (units of length).

`J` = Jacobian matrix
`loc` = location of the quadrature point in physical coordinates,
`conn` = connectivity of the element,
`N` = matrix of basis function values at the quadrature point.
"""
function Jacobiansurface(self::IntegData{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {T<:FESet0Manifold, CC}
    Jac = Jacobianpoint(self, J, loc, conn, N)::FFlt
    if self.axisymmetric
        return Jac*2*pi*loc[1]*self.otherdimension(loc, conn,  N);
    else
        return Jac*self.otherdimension(loc, conn,  N)
    end
end

"""
    Jacobianvolume(self::IntegData{T}, J::FFltMat,
          loc::FFltMat, conn::CC,
          N::FFltMat)::FFlt where {T<:FESet0Manifold, CC}

Evaluate the volume Jacobian.

For the zero-dimensional cell,  the volume Jacobian is
    (i) the product of the point Jacobian and the other dimension
    (units of length cubed);
    or,  when used as axially symmetric
    (ii) the product of the point Jacobian and the circumference of
    the circle through the point `loc` and the other dimension (units
    of length squared).

`J` = Jacobian matrix
`loc` = location of the quadrature point in physical coordinates,
`conn` = connectivity of the element,
`N` = matrix of basis function values at the quadrature point.
"""
function Jacobianvolume(self::IntegData{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {T<:FESet0Manifold, CC}
    Jac = Jacobianpoint(self, J, loc, conn, N)
    if self.axisymmetric
        return Jac*2*pi*loc[1]*self.otherdimension(loc, conn,  N);
    else
        return Jac*self.otherdimension(loc, conn,  N)
    end
end

"""
    Jacobianmdim(self::IntegData{T}, J::FFltMat,
      loc::FFltMat, conn::CC,
      N::FFltMat, m::FInt)::FFlt where {T<:FESet0Manifold, CC}

Evaluate the manifold Jacobian for an m-dimensional manifold.

For an 0-dimensional finite element,  the manifold Jacobian is for
    m=0: +1
    m=1: Jacobiancurve
    m=2: Jacobiansurface
    m=3: Jacobianvolume
"""
function Jacobianmdim(self::IntegData{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat, m::FInt)::FFlt where {T<:FESet0Manifold, CC}
    @assert (m >= 0)  && (m <= 3)
    if (m==3)
        return Jacobianvolume(self, J, loc, conn, N)
    elseif (m==2)
        return Jacobiansurface(self, J, loc, conn, N)
    elseif (m==1)
        return Jacobiancurve(self, J, loc, conn, N)
    else # (m==0)
        return Jacobianpoint(self, J, loc, conn, N)
    end
end


"""
    Jacobiancurve(self::IntegData{T}, J::FFltMat,
              loc::FFltMat, conn::CC,
              N::FFltMat)::FFlt where {T<:FESet1Manifold, CC}

Evaluate the curve Jacobian.

`J` = Jacobian matrix
`loc` = location of the quadrature point in physical coordinates,
`conn` = connectivity of the element,
`N` = matrix of basis function values at the quadrature point.
"""
function Jacobiancurve(self::IntegData{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {T<:FESet1Manifold, CC}
    return Jacobian(self.fes, J)::FFlt
end

"""
    Jacobiansurface(self::IntegData{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                N::FFltMat)::FFlt where {T<:FESet1Manifold, CC}

Evaluate the surface Jacobian.

For the one-dimensional cell,  the surface Jacobian is
    (i) the product of the curve Jacobian and the other dimension
    (units of length);
    or,  when used as axially symmetric
    (ii) the product of the curve Jacobian and the circumference of
    the circle through the point `loc`.

`J` = Jacobian matrix
`loc` = location of the quadrature point in physical coordinates,
`conn` = connectivity of the element,
`N` = matrix of basis function values at the quadrature point.
"""
function Jacobiansurface(self::IntegData{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {T<:FESet1Manifold, CC}
    Jac = Jacobiancurve(self, J, loc, conn, N)
    if self.axisymmetric
        return Jac*2*pi*loc[1];
    else
        return Jac*self.otherdimension(loc, conn,  N)
    end
end

"""
    Jacobianvolume(self::IntegData{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                N::FFltMat)::FFlt where {T<:FESet1Manifold, CC}

Evaluate the volume Jacobian.

For the one-dimensional cell,  the volume Jacobian is
    (i) the product of the curve Jacobian and the other dimension
    (units of length squared);
    or,  when used as axially symmetric
    (ii) the product of the curve Jacobian and the circumference of
    the circle through the point `loc` and the other dimension (units of length).

`J` = Jacobian matrix
`loc` = location of the quadrature point in physical coordinates,
`conn` = connectivity of the element,
`N` = matrix of basis function values at the quadrature point.
"""
function Jacobianvolume(self::IntegData{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {T<:FESet1Manifold, CC}
    Jac = Jacobiancurve(self, J, loc, conn, N)
    if self.axisymmetric
        return Jac*2*pi*loc[1]*self.otherdimension(loc, conn,  N);
    else
        return Jac*self.otherdimension(loc, conn,  N)
    end
end

"""
    Jacobianmdim(self::IntegData{T}, J::FFltMat,
    loc::FFltMat, conn::CC,
    N::FFltMat, m::FInt)::FFlt where {T<:FESet1Manifold, CC}

Evaluate the manifold Jacobian for an m-dimensional manifold.

For an 1-dimensional finite element,  the manifold Jacobian is for
    m=1: Jacobiancurve
    m=2: Jacobiansurface
    m=3: Jacobianvolume
"""
function Jacobianmdim(self::IntegData{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat, m::FInt)::FFlt where {T<:FESet1Manifold, CC}
    @assert (m >= 1) && (m <= 3)
    if (m==3)
        return Jacobianvolume(self, J, loc, conn, N)
    elseif (m==2)
        return Jacobiansurface(self, J, loc, conn, N)
    else # (m==1)
        return Jacobiancurve(self, J, loc, conn, N)
    end
end


"""
    Jacobiansurface(self::IntegData{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                 N::FFltMat)::FFlt where {T<:FESet2Manifold, CC}

Evaluate the surface Jacobian.

`J` = Jacobian matrix
`loc` = location of the quadrature point in physical coordinates,
`conn` = connectivity of the element,
`N` = matrix of basis function values at the quadrature point.
"""
function Jacobiansurface(self::IntegData{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {T<:FESet2Manifold, CC}
    return Jacobian(self.fes, J)
end

"""
    Jacobianvolume(self::IntegData{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                N::FFltMat)::FFlt where {T<:FESet2Manifold, CC}

Evaluate the volume Jacobian.

For the two-dimensional cell,  the volume Jacobian is
    (i) the product of the surface Jacobian and the other dimension
    (units of length);
    or,  when used as axially symmetric
    (ii) the product of the surface Jacobian and the circumference of
    the circle through the point `loc` (units of length).

`J` = Jacobian matrix
`loc` = location of the quadrature point in physical coordinates,
`conn` = connectivity of the element,
`N` = matrix of basis function values at the quadrature point.
"""
function Jacobianvolume(self::IntegData{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {T<:FESet2Manifold, CC}
    Jac = Jacobiansurface(self, J, loc, conn, N)::FFlt
    if self.axisymmetric
        return Jac*2*pi*loc[1];
    else
        return Jac*self.otherdimension(loc, conn,  N)
    end
end

"""
    Jacobianmdim(self::IntegData{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                N::FFltMat, m::FInt)::FFlt where {T<:FESet2Manifold, CC}

Evaluate the manifold Jacobian for an m-dimensional manifold.

For an 2-dimensional finite element,  the manifold Jacobian is for
    m=2: Jacobiansurface
    m=3: Jacobianvolume
"""
function Jacobianmdim(self::IntegData{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat, m::FInt)::FFlt where {T<:FESet2Manifold, CC}
    @assert (m >= 2) && (m <= 3)
    if (m==3)
        return Jacobianvolume(self, J, loc, conn, N)
    else # (m==2)
        return Jacobiansurface(self, J, loc, conn, N)
    end
end


"""
    Jacobianvolume(self::IntegData{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                N::FFltMat)::FFlt where {T<:FESet3Manifold, CC}
Evaluate the volume Jacobian.

`J` = Jacobian matrix
`loc` = location of the quadrature point in physical coordinates,
`conn` = connectivity of the element,
`N` = matrix of basis function values at the quadrature point.
"""
function Jacobianvolume(self::IntegData{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {T<:FESet3Manifold, CC}
    return Jacobian(self.fes, J)::FFlt
end

"""
    Jacobianmdim{T<:FESet3Manifold}(self::IntegData, J::FFltMat,
                loc::FFltMat, conn::FIntMat, N::FFltMat, m::FInt)

Evaluate the manifold Jacobian for an m-dimensional manifold.

For an 3-dimensional cell,  the manifold Jacobian is
    m=3: Jacobianvolume
"""
function Jacobianmdim(self::IntegData{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat, m::FInt)::FFlt where {T<:FESet3Manifold, CC}
    @assert (m == 3)
    return Jacobianvolume(self, J, loc, conn, N)
end

"""
    integrationdata(self::IntegData)

Calculate the data needed for  numerical quadrature.

`npts`, `Ns`, `gradNparams`, `w`, `pc` = number of quadrature points, arrays of
basis function values at the quadrature points,  arrays of gradients of basis
functions  with respect  to the parametric coordinates, array of weights and
array of locations of the quadrature points.
"""
function  integrationdata(self::IntegData)
    pc::FFltMat = self.integration_rule.param_coords;
    w::FFltMat  =  self.integration_rule.weights ;
    npts::FInt = self.integration_rule.npts;
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    Ns = Array{FFltMat}(1,npts);
    gradNparams = Array{FFltMat}(1,npts);
    for j=1:npts
        Ns[j] = bfun(self.fes,vec(pc[j,:]));
        gradNparams[j] = bfundpar(self.fes,vec(pc[j,:]));
    end
    return npts, Ns, gradNparams, w, pc
end


end
