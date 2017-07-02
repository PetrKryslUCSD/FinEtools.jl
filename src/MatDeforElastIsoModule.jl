module MatDeforElastIsoModule

using FinEtools.FTypesModule
using FinEtools.DeforModelRedModule
using FinEtools.MatDeforModule

const mI = diagm([1.0, 1.0, 1.0, 0.5, 0.5, 0.5]);
const m1 = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];

"""
    MatDeforElastIso

Linear isotropic elasticity  material.

```
tangentmoduli!::MTAN
update!::MUPD
thermalstrain!::MTHS
```
"""
struct  MatDeforElastIso{MR<:DeforModelRed,
  MTAN<:Function, MUPD<:Function, MTHS<:Function} <: MatDefor
  mr::Type{MR}
  mass_density::FFlt # mass density
  E::FFlt # Young's modulus
  nu::FFlt # Poisson ratio
  CTE::FFlt # Coefficient of Thermal Expansion
  D::FFltMat # cached matrix of 3D tangent moduli
  tangentmoduli!::MTAN
  update!::MUPD
  thermalstrain!::MTHS
end
export MatDeforElastIso

function threedD(E::FFlt, nu::FFlt)
    lambda = E * nu / (1 + nu) / (1 - 2*(nu));
  mu = E / 2. / (1+nu);
  D = lambda * m1 * m1' + 2. * mu * mI;
  return D
end

function MatDeforElastIso(mr::Type{DeforModelRed3D},
  mass_density::FFlt, E::FFlt, nu::FFlt, CTE::FFlt)
  return MatDeforElastIso(mr, mass_density, E, nu, CTE, threedD(E, nu),
    tangentmoduli3d!, update3d!, thermalstrain3d!)
end

function MatDeforElastIso(mr::Type{DeforModelRed2DStress},
  mass_density::FFlt, E::FFlt, nu::FFlt, CTE::FFlt)
  return MatDeforElastIso(mr, mass_density, E, nu, CTE, threedD(E, nu),
    tangentmoduli2dstrs!, update2dstrs!, thermalstrain2dstrs!)
end

function MatDeforElastIso(mr::Type{DeforModelRed2DStrain},
  mass_density::FFlt, E::FFlt, nu::FFlt, CTE::FFlt)
  return MatDeforElastIso(mr, mass_density, E, nu, CTE, threedD(E, nu),
    tangentmoduli2dstrn!, update2dstrn!, thermalstrain2dstrn!)
end

function MatDeforElastIso(mr::Type{DeforModelRed2DAxisymm},
  mass_density::FFlt, E::FFlt, nu::FFlt, CTE::FFlt)
  return MatDeforElastIso(mr, mass_density, E, nu, CTE, threedD(E, nu),
    tangentmoduli2daxi!, update2daxi!, thermalstrain2daxi!)
end

function MatDeforElastIso(mr::Type{MR}, E::FFlt, nu::FFlt) where {MR}
  mass_density = 1.0
  CTE = 0.0
  return MatDeforElastIso(mr, mass_density, E, nu, CTE)
end

################################################################################
# 3-D solid model
################################################################################

"""
    tangentmoduli3d!(self::MatDeforElastIso,
      ms::MatDeforElastIsoMS, D::FFltMat,
      t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)

Calculate the material stiffness matrix.

`D` = matrix of tangent moduli, 6 x 6, supplied as a buffer and overwritten.
"""
function tangentmoduli3d!(self::MatDeforElastIso,
  D::FFltMat,
  t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)
    copy!(D,self.D);
  return D
end

"""
    update3d!(self::MatDeforElastIso,  stress::FFltVec, output::FFltVec,
      strain::FFltVec, thstrain::FFltVec=zeros(6), t::FFlt= 0.0, dt::FFlt= 0.0,
      loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)

Update material state.

`strain` = strain vector,
`thstrain` = thermal strain vector,
`t` = current time,
`dt` = current time step,
`loc` = location of the quadrature point in global Cartesian coordinates,
`label` = label of the finite element in which the quadrature point is found.

Output:
`stress` = stress vector, allocated by the caller with a size of the number of
  stress and strain components, `nstsstn`. The components of the stress vector are
  calculated and stored in the `stress` vector.
`output` =  array which is (if necessary) allocated  in an appropriate size, filled
  with the output quantity, and returned.
"""
function update3d!(self::MatDeforElastIso,  stress::FFltVec, output::FFltVec,
  strain::FFltVec, thstrain::FFltVec=zeros(6), t::FFlt= 0.0, dt::FFlt= 0.0,
  loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)
  @assert length(stress) == nstsstn(self.mr)
  A_mul_B!(stress, self.D, strain-thstrain);
  if quantity == :nothing
    #Nothing to be copied to the output array
  elseif quantity == :cauchy || quantity == :Cauchy
    (length(output) >= 6) || (output = zeros(6)) # make sure we can store it
    copy!(output, stress);
  elseif quantity == :pressure || quantity == :Pressure
    output[1]  =  -sum(sstress[1:3])/3.
  elseif quantity == :princCauchy || quantity == :princcauchy
    t = zeros(FFlt,3,3)
    t = stress6vto3x3t!(stress,t);
    ep = eig(t);
    (length(output) >= 3) || (output = zeros(3)) # make sure we can store it
    copy!(output[1:3],  sort(ep[1]));
  elseif quantity==:vonMises || quantity==:vonmises || quantity==:von_mises || quantity==:vm
    s1=stress[1]; s2=stress[2]; s3=stress[3];
    s4=stress[4]; s5=stress[5]; s6=stress[6];
    (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
    output[1] = sqrt(1.0/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)))
  end
  return output
end

"""
    thermalstrain3d!(self::MatDeforElastIso, thstrain::FFltMat, dT= 0.0)

Compute thermal strain from the supplied temperature increment.
"""
function thermalstrain3d!(self::MatDeforElastIso, thstrain::FFltVec, dT= 0.0)
  @assert length(thstrain) == nthstn(self.mr)
  thstrain[1] = self.CTE*dT
  thstrain[2] = self.CTE*dT
  thstrain[3] = self.CTE*dT
  thstrain[4] = 0.0
  thstrain[5] = 0.0
  thstrain[6] = 0.0
  return thstrain
end


################################################################################
# 2-D plane stress
################################################################################

"""
    tangentmoduli2dstrs!(self::MatDeforElastIso,
      ms::MatDeforElastIsoMS, D::FFltMat,
      t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)

Calculate the material stiffness matrix.

`D` = matrix of tangent moduli, 3 x 3, supplied as a buffer and overwritten.
"""
function tangentmoduli2dstrs!(self::MatDeforElastIso,
  D::FFltMat,
  t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)
  D[1:2, 1:2] = self.D[1:2, 1:2] -
    (reshape(self.D[1:2,3], 2, 1) * reshape(self.D[3,1:2], 1, 2))/self.D[3, 3]
  const ix=[1, 2, 4];
  for i = 1:3
    D[3,i] = D[i,3] = self.D[4, ix[i]];
  end
  return D
end

"""
    update2dstrs!(self::MatDeforElastIso, output::FFltVec,
      strain::FFltVec, thstrain::FFltVec=zeros(3), t::FFlt= 0.0, dt::FFlt= 0.0,
      loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)

Update material state.

`strain` = strain vector,
`thstrain` = thermal strain vector,
`t` = current time,
`dt` = current time step,
`loc` = location of the quadrature point in global Cartesian coordinates,
`label` = label of the finite element in which the quadrature point is found.

Output:
`stress` = stress vector, allocated by the caller with a size of the number of
  stress and strain components, `nstsstn`. The components of the stress vector are
  calculated and stored in the `stress` vector.
`output` =  array which is (if necessary) allocated  in an appropriate size, filled
  with the output quantity, and returned.
"""
function update2dstrs!(self::MatDeforElastIso, output::FFltVec,
  strain::FFltVec, thstrain::FFltVec=zeros(3), t::FFlt= 0.0, dt::FFlt= 0.0,
  loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)
  @assert length(stress) == nstsstn(self.mr)
  D = zeros(3, 3)
  tangentmoduli2dstrs!(self, D, t, dt, loc, label)
  A_mul_B!(stress, D, strain-thstrain);
  if quantity == :nothing
    #Nothing to be copied to the output array
  elseif quantity == :cauchy || quantity == :Cauchy
    (length(output) >= 3) || (output = zeros(3)) # make sure we can store it
    copy!(output, stress);
  elseif quantity == :pressure || quantity == :Pressure
    (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
    output[1] = -sum(stress[1:2])/2.
  elseif quantity == :princCauchy || quantity == :princcauchy
    t = zeros(FFlt,2,2)
    t = stress3vto2x2t!(stress,t);
    ep = eig(t);
    (length(output) >= 2) || (output = zeros(2)) # make sure we can store it
    output[1:2] =  sort(ep[1])
  elseif quantity==:vonMises || quantity==:vonmises || quantity==:von_mises || quantity==:vm
    s1=stress[1]; s2=stress[2]; s3=0.0;
    s4=stress[3]; s5=0.0; s6=0.0;
    (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
    output[1] = sqrt(1.0/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)))
  end
  return output
end

"""
    thermalstrain2dstrs!(self::MatDeforElastIso, thstrain::FFltVec, dT= 0.0)

Compute thermal strain from the supplied temperature increment.
"""
function thermalstrain2dstrs!(self::MatDeforElastIso, thstrain::FFltVec, dT= 0.0)
  @assert length(thstrain) == nthstn(self.mr)
  thstrain[1] = self.CTE*dT
  thstrain[2] = self.CTE*dT
  thstrain[3] = 0.0
  return thstrain
end


################################################################################
# 2-D plane strain
################################################################################

"""
    tangentmoduli2dstrn!(self::MatDeforElastIso,
      ms::MatDeforElastIsoMS, D::FFltMat,
      t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)

Calculate the material stiffness matrix.

`D` = matrix of tangent moduli, 3 x 3, supplied as a buffer and overwritten.
"""
function tangentmoduli2dstrn!(self::MatDeforElastIso,
  D::FFltMat,
  t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)
  const ix = [1, 2, 4];
  for i = 1:length(ix)
    for j = 1:length(ix)
      D[j,i] = self.D[ix[j], ix[i]];
    end
  end
  return D
end

"""
    update2dstrn!(self::MatDeforElastIso,  stress::FFltVec, output::FFltVec,
      strain::FFltVec, thstrain::FFltVec=zeros(3), t::FFlt= 0.0, dt::FFlt= 0.0,
      loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)

Update material state.

`strain` = strain vector,
`thstrain` = thermal strain vector,
`t` = current time,
`dt` = current time step,
`loc` = location of the quadrature point in global Cartesian coordinates,
`label` = label of the finite element in which the quadrature point is found.

Output:
`stress` = stress vector, allocated by the caller with a size of the number of
  stress and strain components, `nstsstn`. The components of the stress vector are
  calculated and stored in the `stress` vector.
`output` =  array which is (if necessary) allocated  in an appropriate size, filled
  with the output quantity, and returned.
"""
function update2dstrn!(self::MatDeforElastIso,  stress::FFltVec, output::FFltVec,
  strain::FFltVec, thstrain::FFltVec=zeros(3), t::FFlt= 0.0, dt::FFlt= 0.0,
  loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)
  @assert length(stress) == nstsstn(self.mr)
  D = zeros(3, 3)
  tangentmoduli2dstrn!(self, D, t, dt, loc, label)
  A_mul_B!(stress, D, strain-thstrain[1:3]);
  if quantity == :nothing
    #Nothing to be copied to the output array
  elseif quantity == :cauchy || quantity == :Cauchy
    # sigmax, sigmay, tauxy, sigmaz
    # thstrain[4] =The through the thickness thermal strain
    sz = dot(self.D[3, 1:2], strain[1:2]-thstrain[1:2])-self.D[3,3]*thstrain[4];
    (length(output) >= 4) || (output = zeros(4)) # make sure we can store it
    copy!(output, stress); output[4] = sz
  elseif quantity == :pressure || quantity == :Pressure
    (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
    output[1] = -sum(stress[[1,2,4]])/3.
  elseif quantity == :princCauchy || quantity == :princcauchy
    (length(output) >= 3) || (output = zeros(3)) # make sure we can store it
    t = zeros(FFlt, 3,3)
    t = stress4vto3x3t!(t, stress[[1,2,4,3]]);
    ep = eig(t);
    (length(output) >= 3) || (output = zeros(3)) # make sure we can store it
    output[1:3] = sort(ep[1])
  elseif quantity==:vonMises || quantity==:vonmises || quantity==:von_mises || quantity==:vm
    (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
    s1=stress[1]; s2=stress[2]; s3=stress[4];
    s4=stress[3]; s5=0.0; s6=0.0;
    (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
    output[1] = sqrt(1.0/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)))
  end
  return output
end

"""
    thermalstrain2dstrn!(self::MatDeforElastIso, thstrain::FFltVec, dT= 0.0)

Compute thermal strain from the supplied temperature increment.

The thermal strain is evaluated  for the  two normal strains, the shear
strain, and finally for the through the thickness strain.
thstrain = [ex, ey, 0.0, ez].
"""
function thermalstrain2dstrn!(self::MatDeforElastIso, thstrain::FFltVec, dT= 0.0)
  @assert length(thstrain) == nthstn(self.mr)
  thstrain[1] = self.CTE*dT
  thstrain[2] = self.CTE*dT
  thstrain[3] = 0.0
  thstrain[4] = self.CTE*dT
  return thstrain
end

################################################################################
# 2-D axially symmetric
################################################################################

"""
    tangentmoduli2daxi!(self::MatDeforElastIso,
      ms::MatDeforElastIsoMS, D::FFltMat,
      t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)

Calculate the material stiffness matrix.

`D` = matrix of tangent moduli, 3 x 3, supplied as a buffer and overwritten.
"""
function tangentmoduli2daxi!(self::MatDeforElastIso,
  D::FFltMat,
  t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)
  for i = 1:4
    for j = 1:4
      D[j,i] = self.D[i, j];
    end
  end
  return D
end

"""
    update2daxi!(self::MatDeforElastIso,  stress::FFltVec, output::FFltVec,
      strain::FFltVec, thstrain::FFltVec=zeros(3), t::FFlt= 0.0, dt::FFlt= 0.0,
      loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)

Update material state.

`strain` = strain vector,
`thstrain` = thermal strain vector,
`t` = current time,
`dt` = current time step,
`loc` = location of the quadrature point in global Cartesian coordinates,
`label` = label of the finite element in which the quadrature point is found.

Output:
`stress` = stress vector, allocated by the caller with a size of the number of
  stress and strain components, `nstsstn`. The components of the stress vector are
  calculated and stored in the `stress` vector.
`output` =  array which is (if necessary) allocated  in an appropriate size, filled
  with the output quantity, and returned.
"""
function update2daxi!(self::MatDeforElastIso,  stress::FFltVec, output::FFltVec,
  strain::FFltVec, thstrain::FFltVec=zeros(3), t::FFlt= 0.0, dt::FFlt= 0.0,
  loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)
  @assert length(stress) == nstsstn(self.mr)
  D = zeros(4, 4)
  tangentmoduli2daxi!(self, D, t, dt, loc, label)
  A_mul_B!(stress, D, strain-thstrain);
  if quantity == :nothing
    #Nothing to be copied to the output array
  elseif quantity == :cauchy || quantity == :Cauchy
    (length(output) >= 4) || (output = zeros(4)) # make sure we can store it
    copy!(output, stress);
  elseif quantity == :pressure || quantity == :Pressure
    (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
    output[1] = -sum(stress[[1,2,3]])/3.
  elseif quantity == :princCauchy || quantity == :princcauchy
    t = zeros(FFlt,3,3)
    t = stress4vto3x3t!(stress,t);
    ep = eig(t);
    (length(output) >= 3) || (output = zeros(3)) # make sure we can store it
    output[1:3] = sort(ep[1])
  elseif quantity==:vonMises || quantity==:vonmises || quantity==:von_mises || quantity==:vm
    s1=stress[1]; s2=stress[2]; s3=stress[3];
    s4=stress[4]; s5=0.0; s6=0.0;
    (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
    output[1] = sqrt(1.0/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)))
  end
  return output
end

"""
    thermalstrain2daxi!(self::MatDeforElastIso, thstrain::FFltVec, dT= 0.0)

Compute thermal strain from the supplied temperature increment.

The thermal strain is evaluated  for the  three normal strains and the shear
strain.
"""
function thermalstrain2daxi!(self::MatDeforElastIso, thstrain::FFltVec, dT= 0.0)
  @assert length(thstrain) == nthstn(self.mr)
  thstrain[1] = self.CTE*dT
  thstrain[2] = self.CTE*dT
  thstrain[3] = self.CTE*dT
  thstrain[4] = 0.0
  return thstrain
end

end

# ################################################################################
# # 1D model
#
# function tangentmoduli!{P<:PropertyDeformationLinear}(::Type{DeforModelRed1D},
#                         self::MatDeformationLinear{P},
#                         D::FFltMat; context...)
#     # # Calculate the material stiffness matrix.
#     # #
#     # # Arguments
#     # #     m=material
#     # #     context=structure with mandatory and optional fields that are required
#     # # by the specific material implementation.
#     # #
#     # # the output arguments are
#     # #     D=matrix 6x6 in the local material orientation Cartesian basis
#     # #
#
#     D3d=zeros(FFlt,6,6)
#     tangentmoduli3d!(self.property, D3d; context...);
#     D[1,1] = D3d[1, 1]- D3d[1,2:3]*D3d[2:3,2:3]\D3d[2:3,1];
#     return D
# end
# export tangentmoduli!
