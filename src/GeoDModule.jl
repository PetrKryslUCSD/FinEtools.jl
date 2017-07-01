module GeoDModule

using FinEtools
using FinEtools.FTypesModule
using FinEtools.FESetModule
using FinEtools.IntegRuleModule
using FinEtools.NodalFieldModule
using FinEtools.ForceIntensityModule
using FinEtools.AssemblyModule
using FinEtools.CSysModule

"""
    GeoD{T<:FESet}

Geometry data for all finite element models.

T = type of finite element set.  The type of the FE set will be dependent upon the
operations required. For instance, for interior (volume) integrals such as body
load or the stiffness hexahedral H8 may be used whereas for boundary  (surface)
integrals quadrilateral Q4 would be needed.
"""
mutable struct GeoD{S<:FESet, F<:Function}
    fes::S # finite element set object
    integration_rule::IntegRule  # integration rule object
    mcsys::CSys # updater of the material orientation matrix
    otherdimension::F
    axisymmetric::Bool
end
export GeoD

"""
    GeoD(fes::S, integration_rule::IntegRule) where {S<:FESet}

Construct with the default orientation matrix (identity), and the other
dimension  being the default 1.0.
"""
function  GeoD(fes::S, integration_rule::IntegRule) where {S<:FESet}
  return GeoD(fes, integration_rule, CSys(manifdim(fes)),
            otherdimensionunity, false)
end

"""
    GeoD(fes::S, integration_rule::IntegRule,
      otherdimension::FFlt) where {S<:FESet}

Construct with the default orientation matrix (identity), and constant other
dimension.
"""
function  GeoD(fes::S, integration_rule::IntegRule,
  otherdimension::FFlt) where {S<:FESet}
  function otherdimensionfu(loc::FFltMat,
    conn::CC, N::FFltMat)::FFlt where {CC<:AbstractArray{FInt}}
    return otherdimension::FFlt
  end
  return GeoD(fes, integration_rule, CSys(manifdim(fes)),
            otherdimensionfu, false)
end

"""
    GeoD(fes::S, integration_rule::IntegRule,
      axisymmetric::Bool) where {S<:FESet}

Construct with the default orientation matrix (identity), for axially symmetric
models. The other dimension is  the default  unity (1.0).
"""
function GeoD(fes::S, integration_rule::IntegRule,
  axisymmetric::Bool) where {S<:FESet}
  return GeoD(fes, integration_rule, CSys(manifdim(fes)),
            otherdimensionunity, axisymmetric)
end

"""
    GeoD(fes::S, integration_rule::IntegRule,
      mcsys::Csys) where {S<:FESet}

Construct with specified orientation matrix (identity).
"""
function  GeoD(fes::S, integration_rule::IntegRule,
  mcsys::CSys) where {S<:FESet}
  return GeoD(fes, integration_rule, mcsys,
            otherdimensionunity, false)
end

"""
    otherdimensionunity(loc::FFltMat, conn::CC, N::FFltMat)::FFlt
              where {CC<:AbstractArray{FInt}}

Evaluate the other dimension: default is 1.0.
"""
function otherdimensionunity(loc::FFltMat,
  conn::CC, N::FFltMat)::FFlt where {CC<:AbstractArray{FInt}}
  return 1.0
end


"""
    Jacobianpoint(self::GeoD{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                N::FFltMat)::FFlt where {T<:FESet0Manifold, CC<:AbstractArray{FInt}}

Evaluate the point Jacobian.

`J` = Jacobian matrix
`loc` = location of the quadrature point in physical coordinates,
`conn` = connectivity of the element,
`N` = matrix of basis function values at the quadrature point.
"""
function Jacobianpoint(self::GeoD{T}, J::FFltMat,
            loc::FFltMat, conn::CC,
            N::FFltMat)::FFlt where {T<:FESet0Manifold, CC<:AbstractArray{FInt}}
  return FinEtools.FESetModule.Jacobian(self.fes, J)::FFlt
end

"""
    Jacobiancurve{T<:FESet0Manifold}(self::GeoD, J::FFltMat,
                loc::FFltMat, conn::FIntMat, N::FFltMat)

Evaluate the curve Jacobian.

`J` = Jacobian matrix
`loc` = location of the quadrature point in physical coordinates,
`conn` = connectivity of the element,
`N` = matrix of basis function values at the quadrature point.
"""
function Jacobiancurve(self::GeoD{T}, J::FFltMat,
            loc::FFltMat, conn::CC,
            N::FFltMat)::FFlt where {T<:FESet0Manifold, CC<:AbstractArray{FInt}}
  Jac = Jacobianpoint(self, J, loc, conn, N)
  if self.axisymmetric
    return Jac*2*pi*loc[1];
  else
    return Jac*self.otherdimension(loc, conn,  N)
  end
end

"""
    Jacobiansurface(self::GeoD{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                N::FFltMat)::FFlt where {T<:FESet0Manifold, CC<:AbstractArray{FInt}}

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
function Jacobiansurface(self::GeoD{T}, J::FFltMat,
            loc::FFltMat, conn::CC,
            N::FFltMat)::FFlt where {T<:FESet0Manifold, CC<:AbstractArray{FInt}}
  Jac = Jacobianpoint(self, J, loc, conn, N)::FFlt
  if self.axisymmetric
    return Jac*2*pi*loc[1]*self.otherdimension(loc, conn,  N);
  else
    return Jac*self.otherdimension(loc, conn,  N)
  end
end

"""
    Jacobianvolume(self::GeoD{T}, J::FFltMat,
          loc::FFltMat, conn::CC,
          N::FFltMat)::FFlt where {T<:FESet0Manifold, CC<:AbstractArray{FInt}}

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
function Jacobianvolume(self::GeoD{T}, J::FFltMat,
            loc::FFltMat, conn::CC,
            N::FFltMat)::FFlt where {T<:FESet0Manifold, CC<:AbstractArray{FInt}}
  Jac = Jacobianpoint(self, J, loc, conn, N)
  if self.axisymmetric
    return Jac*2*pi*loc[1]*self.otherdimension(loc, conn,  N);
  else
    return Jac*self.otherdimension(loc, conn,  N)
  end
end

"""
    Jacobianmdim(self::GeoD{T}, J::FFltMat,
      loc::FFltMat, conn::CC,
      N::FFltMat, m::FInt)::FFlt where {T<:FESet0Manifold, CC<:AbstractArray{FInt}}

Evaluate the manifold Jacobian for an m-dimensional manifold.

For an 0-dimensional finite element,  the manifold Jacobian is for
    m=0: +1
    m=1: Jacobiancurve
    m=2: Jacobiansurface
    m=3: Jacobianvolume
"""
function Jacobianmdim(self::GeoD{T}, J::FFltMat,
  loc::FFltMat, conn::CC,
  N::FFltMat, m::FInt)::FFlt where {T<:FESet0Manifold, CC<:AbstractArray{FInt}}
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
    Jacobiancurve(self::GeoD{T}, J::FFltMat,
              loc::FFltMat, conn::CC,
              N::FFltMat)::FFlt where {T<:FESet1Manifold, CC<:AbstractArray{FInt}}

Evaluate the curve Jacobian.

`J` = Jacobian matrix
`loc` = location of the quadrature point in physical coordinates,
`conn` = connectivity of the element,
`N` = matrix of basis function values at the quadrature point.
"""
function Jacobiancurve(self::GeoD{T}, J::FFltMat,
          loc::FFltMat, conn::CC,
          N::FFltMat)::FFlt where {T<:FESet1Manifold, CC<:AbstractArray{FInt}}
  return FinEtools.FESetModule.Jacobian(self.fes, J)::FFlt
end

"""
    Jacobiansurface(self::GeoD{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                N::FFltMat)::FFlt where {T<:FESet1Manifold, CC<:AbstractArray{FInt}}

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
function Jacobiansurface(self::GeoD{T}, J::FFltMat,
            loc::FFltMat, conn::CC,
            N::FFltMat)::FFlt where {T<:FESet1Manifold, CC<:AbstractArray{FInt}}
  Jac = Jacobiancurve(self, J, loc, conn, N)
  if self.axisymmetric
    return Jac*2*pi*loc[1];
  else
    return Jac*self.otherdimension(loc, conn,  N)
  end
end

"""
    Jacobianvolume(self::GeoD{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                N::FFltMat)::FFlt where {T<:FESet1Manifold, CC<:AbstractArray{FInt}}

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
function Jacobianvolume(self::GeoD{T}, J::FFltMat,
            loc::FFltMat, conn::CC,
            N::FFltMat)::FFlt where {T<:FESet1Manifold, CC<:AbstractArray{FInt}}
  Jac = Jacobiancurve(self, J, loc, conn, N)
  if self.axisymmetric
    return Jac*2*pi*loc[1]*self.otherdimension(loc, conn,  N);
  else
    return Jac*self.otherdimension(loc, conn,  N)
  end
end

"""
    Jacobianmdim(self::GeoD{T}, J::FFltMat,
    loc::FFltMat, conn::CC,
    N::FFltMat, m::FInt)::FFlt where {T<:FESet1Manifold, CC<:AbstractArray{FInt}}

Evaluate the manifold Jacobian for an m-dimensional manifold.

For an 1-dimensional finite element,  the manifold Jacobian is for
    m=1: Jacobiancurve
    m=2: Jacobiansurface
    m=3: Jacobianvolume
"""
function Jacobianmdim(self::GeoD{T}, J::FFltMat,
      loc::FFltMat, conn::CC,
      N::FFltMat, m::FInt)::FFlt where {T<:FESet1Manifold, CC<:AbstractArray{FInt}}
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
    Jacobiansurface(self::GeoD{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                 N::FFltMat)::FFlt where {T<:FESet2Manifold, CC<:AbstractArray{FInt}}

Evaluate the surface Jacobian.

`J` = Jacobian matrix
`loc` = location of the quadrature point in physical coordinates,
`conn` = connectivity of the element,
`N` = matrix of basis function values at the quadrature point.
"""
function Jacobiansurface(self::GeoD{T}, J::FFltMat,
            loc::FFltMat, conn::CC,
            N::FFltMat)::FFlt where {T<:FESet2Manifold, CC<:AbstractArray{FInt}}
  return FinEtools.FESetModule.Jacobian(self.fes, J)
end

"""
    Jacobianvolume(self::GeoD{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                N::FFltMat)::FFlt where {T<:FESet2Manifold, CC<:AbstractArray{FInt}}

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
function Jacobianvolume(self::GeoD{T}, J::FFltMat,
            loc::FFltMat, conn::CC,
            N::FFltMat)::FFlt where {T<:FESet2Manifold, CC<:AbstractArray{FInt}}
  Jac = Jacobiansurface(self, J, loc, conn, N)::FFlt
  if self.axisymmetric
    return Jac*2*pi*loc[1];
  else
    return Jac*self.otherdimension(loc, conn,  N)
  end
end

"""
    Jacobianmdim(self::GeoD{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                N::FFltMat, m::FInt)::FFlt where {T<:FESet2Manifold, CC<:AbstractArray{FInt}}

Evaluate the manifold Jacobian for an m-dimensional manifold.

For an 2-dimensional finite element,  the manifold Jacobian is for
    m=2: Jacobiansurface
    m=3: Jacobianvolume
"""
function Jacobianmdim(self::GeoD{T}, J::FFltMat,
            loc::FFltMat, conn::CC,
            N::FFltMat, m::FInt)::FFlt where {T<:FESet2Manifold, CC<:AbstractArray{FInt}}
  @assert (m >= 2) && (m <= 3)
  if (m==3)
    return Jacobianvolume(self, J, loc, conn, N)
  else # (m==2)
    return Jacobiansurface(self, J, loc, conn, N)
  end
end


"""
    Jacobianvolume(self::GeoD{T}, J::FFltMat,
                loc::FFltMat, conn::CC,
                N::FFltMat)::FFlt where {T<:FESet3Manifold, CC<:AbstractArray{FInt}}
Evaluate the volume Jacobian.

`J` = Jacobian matrix
`loc` = location of the quadrature point in physical coordinates,
`conn` = connectivity of the element,
`N` = matrix of basis function values at the quadrature point.
"""
function Jacobianvolume(self::GeoD{T}, J::FFltMat,
            loc::FFltMat, conn::CC,
            N::FFltMat)::FFlt where {T<:FESet3Manifold, CC<:AbstractArray{FInt}}
  return FinEtools.FESetModule.Jacobian(self.fes, J)::FFlt
end

"""
    Jacobianmdim{T<:FESet3Manifold}(self::GeoD, J::FFltMat,
                loc::FFltMat, conn::FIntMat, N::FFltMat, m::FInt)

Evaluate the manifold Jacobian for an m-dimensional manifold.

For an 3-dimensional cell,  the manifold Jacobian is
    m=3: Jacobianvolume
"""
function Jacobianmdim(self::GeoD{T}, J::FFltMat,
            loc::FFltMat, conn::CC,
            N::FFltMat, m::FInt)::FFlt where {T<:FESet3Manifold, CC<:AbstractArray{FInt}}
  @assert (m == 3)
  return Jacobianvolume(self, J, loc, conn, N)
end

export Jacobianpoint
export Jacobiancurve
export Jacobiansurface
export Jacobianvolume
export Jacobianmdim

"""
    integrationdata(self::GeoD)

Calculate the data needed for  numerical quadrature.

`npts`, `Ns`, `gradNparams`, `w`, `pc` = number of quadrature points, arrays of
basis function values at the quadrature points,  arrays of gradients of basis
functions  with respect  to the parametric coordinates, array of weights and
array of locations of the quadrature points.
"""
function  integrationdata(self::GeoD)
    # Calculate the data needed for  numerical quadrature.

    pc::FFltMat = self.integration_rule.param_coords;
    w::FFltMat  =  self.integration_rule.weights ;
    npts::FInt = self.integration_rule.npts;
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    Ns = Array{FFltMat}(1,npts);
    gradNparams = Array{FFltMat}(1,npts);
    for j=1:npts
        Ns[j] = FESetModule.bfun(self.fes,vec(pc[j,:]));
        gradNparams[j] = FESetModule.bfundpar(self.fes,vec(pc[j,:]));
    end
    return npts, Ns, gradNparams, w, pc
end
export integrationdata


#  Associated geometry field. Default is there is none, so any
#  operation that requires a geometry field needs to be supplied it.
#  There may be operations that could benefit from pre-computations
#  that involve a geometry field. If so, associating the geometry
#  field gives the FEMM a chance to save on repeated computations.
#         assoc_geom = []; % associated geometry field
#     end

end
