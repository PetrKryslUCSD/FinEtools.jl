module MatDeforLinearElasticModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.DeforModelRedModule: AbstractDeforModelRed
import FinEtools.MatDeforModule: AbstractMatDefor

"""
    AbstractMatDeforLinearElastic <: AbstractMatDefor

Abstract Linear Elasticity  material.

"""
abstract type AbstractMatDeforLinearElastic <: AbstractMatDefor; end

"""
    tangentmoduli!(self::AbstractMatDeforLinearElastic,  D::FFltMat,  t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)

Calculate the material stiffness matrix.

- `D` = matrix of tangent moduli, supplied as a buffer and overwritten. Returned
as output.
"""
function tangentmoduli!(self::AbstractMatDeforLinearElastic,  D::FFltMat,  t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)
    return self.tangentmoduli!(self, D, t, dt, loc, label)
end

"""
    update!(self::AbstractMatDeforLinearElastic,  stress::FFltVec, output::FFltVec,  strain::FFltVec, thstrain::FFltVec=zeros(6), t::FFlt= 0.0, dt::FFlt= 0.0,  loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)

Update material state.

- `strain` = strain vector,
- `thstrain` = thermal strain vector,
- `t` = current time,
- `dt` = current time step,
- `loc` = location of the quadrature point in global Cartesian coordinates,
- `label` = label of the finite element in which the quadrature point is found.

# Output
- `stress` = stress vector, allocated by the caller with a size of the number of stress and
strain components, `nstressstrain`. The components of the stress vector are
calculated and stored in the `stress` vector.
- `output` =  array which is (if necessary) allocated  in an appropriate size, filled
  with the output quantity, and returned.
"""
function update!(self::AbstractMatDeforLinearElastic,  stress::FFltVec, output::FFltVec,  strain::FFltVec, thstrain::FFltVec=zeros(6), t::FFlt= 0.0, dt::FFlt= 0.0,  loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)
    return self.update!(self, stress, output, strain, thstrain, t, dt, loc, label, quantity)
end

"""
    thermalstrain!(self::AbstractMatDeforLinearElastic, thstrain::FFltVec, dT= 0.0)

Compute thermal strain from the supplied temperature increment.

- `thstrain` = thermal strain vector, supplied as buffer, returned as output. 
"""
function thermalstrain!(self::AbstractMatDeforLinearElastic, thstrain::FFltVec, dT= 0.0)
    return self.thermalstrain!(self, thstrain, dT)
end

end
