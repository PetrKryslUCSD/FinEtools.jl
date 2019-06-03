"""
    IntegDomainModule

Module to manage integration domains.

An integration domain consists of the finite elements that approximate the
geometry, the function to supply the "missing" (other) dimension, indication
whether or not the integration domain represents an axially symmetric
situation, and integration rule used to evaluate integrals over the domain.
"""
module IntegDomainModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FESetModule: AbstractFESet, AbstractFESet0Manifold, AbstractFESet1Manifold, AbstractFESet2Manifold, AbstractFESet3Manifold, Jacobian, bfun, bfundpar
import FinEtools.FENodeSetModule: FENodeSet
import FinEtools.IntegRuleModule: AbstractIntegRule

"""
    IntegDomain{S<:AbstractFESet, F<:Function}

Integration domain.

- `T` = type of finite element set.  The type of the FE set will be dependent
  upon the operations required. For instance, for interior (volume) integrals
  such as body load or the stiffness hexahedral H8 may be used, whereas for
  boundary  (surface) integrals quadrilateral Q4 would be needed.
- `F` = type of function to return the "other" dimension.

An integration domain consists of the finite elements that approximate the
geometry, the function to supply the "missing" (other) dimension, indication
whether or not the integration domain represents an axially symmetric
situation, and integration rule used to evaluate integrals over the domain.
"""
mutable struct IntegDomain{S<:AbstractFESet, F<:Function}
    fes::S # finite element set object
    integration_rule::AbstractIntegRule  # integration rule object
    otherdimension::F # function to compute the "other" dimension (thickness, or cross-sectional area)
    axisymmetric::Bool
end

"""
    IntegDomain(fes::S, integration_rule::AbstractIntegRule) where {S<:AbstractFESet}

Construct with the default orientation matrix (identity), and the other
dimension  being the default 1.0.
"""
function  IntegDomain(fes::S, integration_rule::AbstractIntegRule) where {S<:AbstractFESet}
    return IntegDomain(fes, integration_rule, otherdimensionunity, false)
end

"""
    IntegDomain(fes::S, integration_rule::AbstractIntegRule,
      otherdimension::FFlt) where {S<:AbstractFESet}

Construct with the default orientation matrix (identity), and constant other
dimension.
"""
function  IntegDomain(fes::S, integration_rule::AbstractIntegRule, otherdimension::FFlt) where {S<:AbstractFESet}
    function otherdimensionfu(loc::FFltMat,
        conn::CC, N::FFltMat)::FFlt where {CC}
        return otherdimension::FFlt
    end
    return IntegDomain(fes, integration_rule, otherdimensionfu, false)
end

"""
    IntegDomain(fes::S, integration_rule::AbstractIntegRule,
      axisymmetric::Bool) where {S<:AbstractFESet}

Construct with the default orientation matrix (identity), for axially
symmetric models. The other dimension is the default unity (1.0).

This will probably be called when `axisymmetric = true`, since the default is
`axisymmetric = false`.
"""
function IntegDomain(fes::S, integration_rule::AbstractIntegRule, axisymmetric::Bool) where {S<:AbstractFESet}
    return IntegDomain(fes, integration_rule, otherdimensionunity, axisymmetric)
end

"""
    IntegDomain(fes::S, integration_rule::AbstractIntegRule, axisymmetric::Bool,
      otherdimension::FFlt) where {S<:AbstractFESet}

Construct for axially symmetric models. The other dimension is given as a number.
"""
function IntegDomain(fes::S, integration_rule::AbstractIntegRule, axisymmetric::Bool, otherdimension::FFlt) where {S<:AbstractFESet}
    function otherdimensionfu(loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {CC}
        return otherdimension::FFlt
    end
    return IntegDomain(fes, integration_rule, otherdimensionfu, axisymmetric)
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
    Jacobianpoint(self::IntegDomain{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                N::FFltMat)::FFlt where {T<:AbstractFESet0Manifold, CC}

Evaluate the point Jacobian.

- `J` = Jacobian matrix
- `loc` = location of the quadrature point in physical coordinates,
- `conn` = connectivity of the element,
- `N` = matrix of basis function values at the quadrature point.
"""
function Jacobianpoint(self::IntegDomain{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {T<:AbstractFESet0Manifold, CC}
    return Jacobian(self.fes, J)::FFlt
end

"""
    Jacobiancurve(self::IntegDomain{T}, J::FFltMat,
        loc::FFltMat, conn::CC,
        N::FFltMat)::FFlt where {T<:AbstractFESet0Manifold, CC}

Evaluate the curve Jacobian.

- `J` = Jacobian matrix
- `loc` = location of the quadrature point in physical coordinates,
- `conn` = connectivity of the element,
- `N` = matrix of basis function values at the quadrature point.
"""
function Jacobiancurve(self::IntegDomain{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {T<:AbstractFESet0Manifold, CC}
    Jac = Jacobianpoint(self, J, loc, conn, N)
    if self.axisymmetric
        return Jac*2*pi*loc[1];
    else
        return Jac*self.otherdimension(loc, conn,  N)
    end
end

"""
    Jacobiansurface(self::IntegDomain{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                N::FFltMat)::FFlt where {T<:AbstractFESet0Manifold, CC}

Evaluate the surface Jacobian.

For the zero-dimensional cell, the surface Jacobian is (i) the product of the
point Jacobian and the other dimension (units of length squared); or,  when
used as axially symmetric (ii) the product of the point Jacobian and the
circumference of the circle through the point `loc` times the other dimension
(units of length).

- `J` = Jacobian matrix
- `loc` = location of the quadrature point in physical coordinates,
- `conn` = connectivity of the element,
- `N` = matrix of basis function values at the quadrature point.
"""
function Jacobiansurface(self::IntegDomain{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {T<:AbstractFESet0Manifold, CC}
    Jac = Jacobianpoint(self, J, loc, conn, N)::FFlt
    if self.axisymmetric
        return Jac*2*pi*loc[1]*self.otherdimension(loc, conn,  N);
    else
        return Jac*self.otherdimension(loc, conn,  N)
    end
end

"""
    Jacobianvolume(self::IntegDomain{T}, J::FFltMat,
          loc::FFltMat, conn::CC,
          N::FFltMat)::FFlt where {T<:AbstractFESet0Manifold, CC}

Evaluate the volume Jacobian.

For the zero-dimensional cell, the volume Jacobian is (i) the product of the
point Jacobian and the other dimension (units of length cubed); or,  when used
as axially symmetric (ii) the product of the point Jacobian and the
circumference of the circle through the point `loc` and the other dimension
(units of length squared).

- `J` = Jacobian matrix
- `loc` = location of the quadrature point in physical coordinates,
- `conn` = connectivity of the element,
- `N` = matrix of basis function values at the quadrature point.
"""
function Jacobianvolume(self::IntegDomain{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {T<:AbstractFESet0Manifold, CC}
    Jac = Jacobianpoint(self, J, loc, conn, N)
    if self.axisymmetric
        return Jac*2*pi*loc[1]*self.otherdimension(loc, conn,  N);
    else
        return Jac*self.otherdimension(loc, conn,  N)
    end
end

"""
    Jacobianmdim(self::IntegDomain{T}, J::FFltMat,
      loc::FFltMat, conn::CC,
      N::FFltMat, m::FInt)::FFlt where {T<:AbstractFESet0Manifold, CC}

Evaluate the manifold Jacobian for an m-dimensional manifold.

For an 0-dimensional finite element,  the manifold Jacobian is for
- m=0: +1
- m=1: `Jacobiancurve`
- m=2: `Jacobiansurface`
- m=3: `Jacobianvolume`
"""
function Jacobianmdim(self::IntegDomain{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat, m::FInt)::FFlt where {T<:AbstractFESet0Manifold, CC}
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
    Jacobiancurve(self::IntegDomain{T}, J::FFltMat,
              loc::FFltMat, conn::CC,
              N::FFltMat)::FFlt where {T<:AbstractFESet1Manifold, CC}

Evaluate the curve Jacobian.

- `J` = Jacobian matrix
- `loc` = location of the quadrature point in physical coordinates,
- `conn` = connectivity of the element,
- `N` = matrix of basis function values at the quadrature point.
"""
function Jacobiancurve(self::IntegDomain{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {T<:AbstractFESet1Manifold, CC}
    return Jacobian(self.fes, J)::FFlt
end

"""
    Jacobiansurface(self::IntegDomain{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                N::FFltMat)::FFlt where {T<:AbstractFESet1Manifold, CC}

Evaluate the surface Jacobian.

For the one-dimensional cell,  the surface Jacobian is (i) the product of the
curve Jacobian and the other dimension (units of length); or,  when used as
axially symmetric (ii) the product of the curve Jacobian and the circumference
of the circle through the point `loc`.

- `J` = Jacobian matrix
- `loc` = location of the quadrature point in physical coordinates,
- `conn` = connectivity of the element,
- `N` = matrix of basis function values at the quadrature point.
"""
function Jacobiansurface(self::IntegDomain{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {T<:AbstractFESet1Manifold, CC}
    Jac = Jacobiancurve(self, J, loc, conn, N)
    if self.axisymmetric
        return Jac*2*pi*loc[1];
    else
        return Jac*self.otherdimension(loc, conn,  N)
    end
end

"""
    Jacobianvolume(self::IntegDomain{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                N::FFltMat)::FFlt where {T<:AbstractFESet1Manifold, CC}

Evaluate the volume Jacobian.

For the one-dimensional cell,  the volume Jacobian is (i) the product of the
curve Jacobian and the other dimension (units of length squared); or,  when
used as axially symmetric (ii) the product of the curve Jacobian and the
circumference of the circle through the point `loc` and the other dimension
(units of length).

- `J` = Jacobian matrix
- `loc` = location of the quadrature point in physical coordinates,
- `conn` = connectivity of the element,
- `N` = matrix of basis function values at the quadrature point.
"""
function Jacobianvolume(self::IntegDomain{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {T<:AbstractFESet1Manifold, CC}
    Jac = Jacobiancurve(self, J, loc, conn, N)
    if self.axisymmetric
        return Jac*2*pi*loc[1]*self.otherdimension(loc, conn,  N);
    else
        return Jac*self.otherdimension(loc, conn,  N)
    end
end

"""
    Jacobianmdim(self::IntegDomain{T}, J::FFltMat,
    loc::FFltMat, conn::CC,
    N::FFltMat, m::FInt)::FFlt where {T<:AbstractFESet1Manifold, CC}

Evaluate the manifold Jacobian for an m-dimensional manifold.

For an 1-dimensional finite element,  the manifold Jacobian is for
- m=1: `Jacobiancurve`
- m=2: `Jacobiansurface`
- m=3: `Jacobianvolume`
"""
function Jacobianmdim(self::IntegDomain{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat, m::FInt)::FFlt where {T<:AbstractFESet1Manifold, CC}
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
    Jacobiansurface(self::IntegDomain{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                 N::FFltMat)::FFlt where {T<:AbstractFESet2Manifold, CC}

Evaluate the surface Jacobian.

- `J` = Jacobian matrix
- `loc` = location of the quadrature point in physical coordinates,
- `conn` = connectivity of the element,
- `N` = matrix of basis function values at the quadrature point.
"""
function Jacobiansurface(self::IntegDomain{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {T<:AbstractFESet2Manifold, CC}
    return Jacobian(self.fes, J)
end

"""
    Jacobianvolume(self::IntegDomain{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                N::FFltMat)::FFlt where {T<:AbstractFESet2Manifold, CC}

Evaluate the volume Jacobian.

For the two-dimensional cell,  the volume Jacobian is (i) the product of the
surface Jacobian and the other dimension (units of length); or,  when used as
axially symmetric (ii) the product of the surface Jacobian and the
circumference of the circle through the point `loc` (units of length).

- `J` = Jacobian matrix
- `loc` = location of the quadrature point in physical coordinates,
- `conn` = connectivity of the element,
- `N` = matrix of basis function values at the quadrature point.
"""
function Jacobianvolume(self::IntegDomain{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {T<:AbstractFESet2Manifold, CC}
    Jac = Jacobiansurface(self, J, loc, conn, N)::FFlt
    if self.axisymmetric
        return Jac*2*pi*loc[1];
    else
        return Jac*self.otherdimension(loc, conn,  N)
    end
end

"""
    Jacobianmdim(self::IntegDomain{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                N::FFltMat, m::FInt)::FFlt where {T<:AbstractFESet2Manifold, CC}

Evaluate the manifold Jacobian for an m-dimensional manifold.

For an 2-dimensional finite element,  the manifold Jacobian is for
- m=2: `Jacobiansurface`
- m=3: `Jacobianvolume`
"""
function Jacobianmdim(self::IntegDomain{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat, m::FInt)::FFlt where {T<:AbstractFESet2Manifold, CC}
    @assert (m >= 2) && (m <= 3)
    if (m==3)
        return Jacobianvolume(self, J, loc, conn, N)
    else # (m==2)
        return Jacobiansurface(self, J, loc, conn, N)
    end
end


"""
    Jacobianvolume(self::IntegDomain{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                N::FFltMat)::FFlt where {T<:AbstractFESet3Manifold, CC}
Evaluate the volume Jacobian.

- `J` = Jacobian matrix
- `loc` = location of the quadrature point in physical coordinates,
- `conn` = connectivity of the element,
- `N` = matrix of basis function values at the quadrature point.
"""
function Jacobianvolume(self::IntegDomain{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {T<:AbstractFESet3Manifold, CC}
    return Jacobian(self.fes, J)::FFlt
end

"""
    Jacobianmdim{T<:AbstractFESet3Manifold}(self::IntegDomain, J::FFltMat,
                loc::FFltMat, conn::FIntMat, N::FFltMat, m::FInt)

Evaluate the manifold Jacobian for an m-dimensional manifold.

For an 3-dimensional cell,  the manifold Jacobian is
- m=3: `Jacobianvolume`
"""
function Jacobianmdim(self::IntegDomain{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat, m::FInt)::FFlt where {T<:AbstractFESet3Manifold, CC}
    @assert (m == 3)
    return Jacobianvolume(self, J, loc, conn, N)
end

"""
    integrationdata(self::IntegDomain)

Calculate the data needed for  numerical quadrature for the integration rule
stored by the integration domain.
"""
function  integrationdata(self::IntegDomain)
    return integrationdata(self, self.integration_rule)
end


"""
    integrationdata(self::IntegDomain, integration_rule::T) where {T<:AbstractIntegRule}

Calculate the data needed for numerical quadrature.

For given integration domain, compute the quantities needed for numerical
integration.

# Return
`npts`, `Ns`, `gradNparams`, `w`, `pc` = number of quadrature points, arrays of
basis function values at the quadrature points,  arrays of gradients of basis
functions  with respect  to the parametric coordinates, array of weights and
array of locations of the quadrature points.
"""
function  integrationdata(self::IntegDomain, integration_rule::T) where {T<:AbstractIntegRule}
    pc::FFltMat = integration_rule.param_coords;
    w::FFltMat  =  integration_rule.weights ;
    npts::FInt = integration_rule.npts;
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    Ns = FFltMat[];
    gradNparams = FFltMat[];
    for j=1:npts
        push!(Ns, bfun(self.fes,vec(pc[j,:])));
        push!(gradNparams, bfundpar(self.fes,vec(pc[j,:])));
    end
    return npts, reshape(Ns, 1,npts), reshape(gradNparams, 1,npts), w, pc
end

end
