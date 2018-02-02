"""
    MatDeforElastOrthoModule

Module for  orthotropic elastic material.
"""
module MatDeforElastOrthoModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.DeforModelRedModule: DeforModelRed, DeforModelRed3D, DeforModelRed2DStrain, DeforModelRed2DStress, DeforModelRed2DAxisymm, DeforModelRed1D, nstressstrain, nthermstrain
import FinEtools.MatDeforModule: MatDefor, stress6vto3x3t!, stress3vto2x2t!, stress4vto3x3t!, stress4vto3x3t!
if VERSION < v"0.7-"
    copyto! = copy!
end
if VERSION >= v"0.7-"
    At_mul_B!(C, A, B) = Base.LinAlg.mul!(C, Transpose(A), B)
    A_mul_B!(C, A, B) = Base.LinAlg.mul!(C, A, B)
end

"""
    MatDeforElastOrtho

Linear orthotropic elasticity  material.


```tangentmoduli!(self::MatDeforElastOrtho,
  ms::MatDeforElastOrthoMS, D::FFltMat,
  t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)
```


"""
struct  MatDeforElastOrtho{MR<:DeforModelRed,  MTAN<:Function, MUPD<:Function, MTHS<:Function} <: MatDefor
    mr::Type{MR}
    mass_density::FFlt # mass density
    E1::FFlt # Young's modulus for material direction 1
    E2::FFlt # Young's modulus for material direction 2
    E3::FFlt # Young's modulus for material direction 3
    nu12::FFlt
    nu13::FFlt
    nu23::FFlt
    G12::FFlt
    G13::FFlt
    G23::FFlt
    CTE1::FFlt # three thermal expansion coefficients
    CTE2::FFlt # three thermal expansion coefficients
    CTE3::FFlt # three thermal expansion coefficients
    D::Array{FFlt, 2} # cached matrix of tangent moduli
    tangentmoduli!::MTAN
    update!::MUPD
    thermalstrain!::MTHS
end

function threedD(E1::FFlt, E2::FFlt, E3::FFlt,
    nu12::FFlt, nu13::FFlt, nu23::FFlt,
    G12::FFlt, G13::FFlt, G23::FFlt)
    C = [  1.0/E1      -nu12/E1    -nu13/E1  0.0   0.0   0.0;
        -nu12/E1      1.0/E2      -nu23/E2  0.0   0.0   0.0;
        -nu13/E1    -nu23/E2       1.0/E3   0.0   0.0   0.0;
        0.0           0.0           0.0   1/G12   0.0   0.0;
        0.0           0.0           0.0     0.0 1/G13   0.0;
        0.0           0.0           0.0     0.0   0.0  1/G23];
    D = inv(C)
    if (rank(D)<6)
        error("Non-positive definite D!");
    end
    ev = eigvals(D)
    for  e in ev
        #println("$e")
        if (e < 0.0)
        error("Indefinite D!");
        end
    end
    return D
end


function MatDeforElastOrtho(mr::Type{DeforModelRed3D},
    mass_density::FFlt, E1::FFlt, E2::FFlt, E3::FFlt,
    nu12::FFlt, nu13::FFlt, nu23::FFlt,
    G12::FFlt, G13::FFlt, G23::FFlt, CTE1::FFlt, CTE2::FFlt, CTE3::FFlt)
    return MatDeforElastOrtho(mr, mass_density, E1, E2, E3,    nu12, nu13, nu23,
        G12, G13, G23, CTE1, CTE2, CTE3,
        threedD(E1, E2, E3,    nu12, nu13, nu23,    G12, G13, G23),
    tangentmoduli3d!, update3d!, thermalstrain3d!)
end

function MatDeforElastOrtho(mr::Type{DeforModelRed2DStress},
    mass_density::FFlt, E1::FFlt, E2::FFlt, E3::FFlt,
    nu12::FFlt, nu13::FFlt, nu23::FFlt,
    G12::FFlt, G13::FFlt, G23::FFlt, CTE1::FFlt, CTE2::FFlt, CTE3::FFlt)
    return MatDeforElastOrtho(mr, mass_density, E1, E2, E3,    nu12, nu13, nu23,
        G12, G13, G23, CTE1, CTE2, CTE3,
        threedD(E1, E2, E3,    nu12, nu13, nu23,    G12, G13, G23),
        tangentmoduli2dstrs!, update2dstrs!, thermalstrain2dstrs!)
end

function MatDeforElastOrtho(mr::Type{DeforModelRed2DStrain},
    mass_density::FFlt, E1::FFlt, E2::FFlt, E3::FFlt,
    nu12::FFlt, nu13::FFlt, nu23::FFlt,
    G12::FFlt, G13::FFlt, G23::FFlt, CTE1::FFlt, CTE2::FFlt, CTE3::FFlt)
    return MatDeforElastOrtho(mr, mass_density, E1, E2, E3,    nu12, nu13, nu23,
        G12, G13, G23, CTE1, CTE2, CTE3,
        threedD(E1, E2, E3,    nu12, nu13, nu23,    G12, G13, G23),
        tangentmoduli2dstrn!, update2dstrn!, thermalstrain2dstrn!)
end

function MatDeforElastOrtho(mr::Type{DeforModelRed2DAxisymm},
    mass_density::FFlt, E1::FFlt, E2::FFlt, E3::FFlt,
    nu12::FFlt, nu13::FFlt, nu23::FFlt,
    G12::FFlt, G13::FFlt, G23::FFlt, CTE1::FFlt, CTE2::FFlt, CTE3::FFlt)
    return MatDeforElastOrtho(mr, mass_density, E1, E2, E3,    nu12, nu13, nu23,
        G12, G13, G23, CTE1, CTE2, CTE3,
        threedD(E1, E2, E3,    nu12, nu13, nu23,    G12, G13, G23),
        tangentmoduli2daxi!, update2daxi!, thermalstrain2daxi!)
end

function MatDeforElastOrtho(mr::Type{MR}, E1::FFlt, E2::FFlt, E3::FFlt,
    nu12::FFlt, nu13::FFlt, nu23::FFlt,
    G12::FFlt, G13::FFlt, G23::FFlt) where {MR}
    mass_density = 1.0
    CTE1 = CTE2 = CTE3 = 0.0
    return MatDeforElastOrtho(mr, mass_density, E1, E2, E3,    nu12, nu13, nu23,
    G12, G13, G23, CTE1, CTE2, CTE3)
end

function MatDeforElastOrtho(mr::Type{MR}, E::FFlt, nu::FFlt) where {MR}
    mass_density = 1.0
    CTE1 = CTE2 = CTE3 = 0.0
    G = E / 2.0 / (1 + nu)
    return MatDeforElastOrtho(mr, mass_density, E, E, E,    nu, nu, nu,    G, G, G, CTE1, CTE2, CTE3)
end

function MatDeforElastOrtho(mr::Type{MR}, mass_density::FFlt,  E::FFlt, nu::FFlt, CTE::FFlt) where {MR}
    CTE1 = CTE2 = CTE3 = CTE
    G = E / 2.0 / (1 + nu)
    return MatDeforElastOrtho(mr, mass_density, E, E, E,    nu, nu, nu,    G, G, G, CTE1, CTE2, CTE3)
end


################################################################################
# 3-D solid model
################################################################################

"""
    tangentmoduli3d!(self::MatDeforElastOrtho,
      ms::MatDeforElastOrthoMS, D::FFltMat,
      t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)

Calculate the material stiffness matrix.

`D` = matrix of tangent moduli, 6 x 6, supplied as a buffer and overwritten.
"""
function tangentmoduli3d!(self::MatDeforElastOrtho,  D::FFltMat,  t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)
    copyto!(D,self.D);
    return D
end

"""
    update3d!(self::MatDeforElastOrtho,  stress::FFltVec, output::FFltVec,
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
  stress and strain components, `nstressstrain`. The components of the stress vector are
  calculated and stored in the `stress` vector.
`output` =  array which is (if necessary) allocated  in an appropriate size, filled
  with the output quantity, and returned.
"""
function update3d!(self::MatDeforElastOrtho,  stress::FFltVec, output::FFltVec,  strain::FFltVec, thstrain::FFltVec=zeros(6), t::FFlt= 0.0, dt::FFlt= 0.0,  loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)
    @assert length(stress) == nstressstrain(self.mr)
    A_mul_B!(stress, self.D, strain-thstrain);
    if quantity == :nothing
        #Nothing to be copied to the output array
    elseif quantity == :cauchy || quantity == :Cauchy
        (length(output) >= 6) || (output = zeros(6)) # make sure we can store it
        copyto!(output, stress);
    elseif quantity == :pressure || quantity == :Pressure
        output[1]  =  -sum(stress[1:3])/3.
    elseif quantity == :princCauchy || quantity == :princcauchy
        t = zeros(FFlt,3,3)
        t = stress6vto3x3t!(t,stress);
        ep = eig(t);
        (length(output) >= 3) || (output = zeros(3)) # make sure we can store it
        copyto!(output,  sort(ep[1], rev=true));
    elseif quantity==:vonMises || quantity==:vonmises || quantity==:von_mises || quantity==:vm
        s1=stress[1]; s2=stress[2]; s3=stress[3];
        s4=stress[4]; s5=stress[5]; s6=stress[6];
        (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
        output[1] = sqrt(1.0/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)))
    end
    return output
end

"""
    thermalstrain3d!(self::MatDeforElastOrtho, thstrain::FFltMat, dT= 0.0)

Compute thermal strain from the supplied temperature increment.
"""
function thermalstrain3d!(self::MatDeforElastOrtho, thstrain::FFltVec, dT= 0.0)
    thstrain[1] = self.CTE1*dT
    thstrain[2] = self.CTE2*dT
    thstrain[3] = self.CTE3*dT
    thstrain[4] = 0.0
    thstrain[5] = 0.0
    thstrain[6] = 0.0
    return thstrain
end

################################################################################
# 2-D plane stress
################################################################################

"""
    tangentmoduli2dstrs!(self::MatDeforElastOrtho,
      D::FFltMat,
      t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)

Calculate the material stiffness matrix.

`D` = matrix of tangent moduli, 3 x 3, supplied as a buffer and overwritten.
"""
function tangentmoduli2dstrs!(self::MatDeforElastOrtho,  D::FFltMat,  t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)
    D[1:2, 1:2] = self.D[1:2, 1:2] - (reshape(self.D[1:2,3], 2, 1) * reshape(self.D[3,1:2], 1, 2))/self.D[3, 3]
    ix=[1, 2, 4];
    for i = 1:3
        D[3,i] = D[i,3] = self.D[4, ix[i]];
    end
    return D
end

"""
    update2dstrs!(self::MatDeforElastOrtho, stress::FFltVec,  output::FFltVec,
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
  stress and strain components, `nstressstrain`. The components of the stress vector are
  calculated and stored in the `stress` vector.
`output` =  array which is (if necessary) allocated  in an appropriate size, filled
  with the output quantity, and returned.
"""
function update2dstrs!(self::MatDeforElastOrtho, stress::FFltVec,  output::FFltVec,  strain::FFltVec, thstrain::FFltVec=zeros(3), t::FFlt= 0.0, dt::FFlt= 0.0,  loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)
    @assert length(stress) == nstressstrain(self.mr)
    D = zeros(3, 3)
    tangentmoduli2dstrs!(self, D, t, dt, loc, label)
    A_mul_B!(stress, D, strain-thstrain);
    if quantity == :nothing
        #Nothing to be copied to the output array
    elseif quantity == :cauchy || quantity == :Cauchy
        (length(output) >= 3) || (output = zeros(3)) # make sure we can store it
        copyto!(output, stress);
    elseif quantity == :pressure || quantity == :Pressure
        (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
        output[1] = -sum(stress[1:2])/3.
    elseif quantity == :princCauchy || quantity == :princcauchy
        t = zeros(FFlt,2,2)
        t = stress3vto2x2t!(t,stress);
        ep = eig(t);
        (length(output) >= 2) || (output = zeros(2)) # make sure we can store it
        copyto!(output,  sort(ep[1], rev=true));
    elseif quantity==:vonMises || quantity==:vonmises || quantity==:von_mises || quantity==:vm
        s1=stress[1]; s2=stress[2]; s3=0.0;
        s4=stress[3]; s5=0.0; s6=0.0;
        (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
        output[1] = sqrt(1.0/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)))
    end
    return output
end

"""
    thermalstrain2dstrs!(self::MatDeforElastOrtho, thstrain::FFltVec, dT= 0.0)

Compute thermal strain from the supplied temperature increment.
"""
function thermalstrain2dstrs!(self::MatDeforElastOrtho, thstrain::FFltVec, dT= 0.0)
    @assert length(thstrain) == nthermstrain(self.mr)
    thstrain[1] = self.CTE1*dT
    thstrain[2] = self.CTE2*dT
    thstrain[3] = 0.0
    return thstrain
end


################################################################################
# 2-D plane strain
################################################################################

"""
    tangentmoduli2dstrn!(self::MatDeforElastOrtho,
      D::FFltMat,
      t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)

Calculate the material stiffness matrix.

`D` = matrix of tangent moduli, 3 x 3, supplied as a buffer and overwritten.
"""
function tangentmoduli2dstrn!(self::MatDeforElastOrtho, D::FFltMat,  t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)
    ix = [1, 2, 4];
    for i = 1:length(ix)
        for j = 1:length(ix)
        D[j,i] = self.D[ix[j], ix[i]];
        end
    end
    return D
end

"""
    update2dstrn!(self::MatDeforElastOrtho,  stress::FFltVec, output::FFltVec,
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
  stress and strain components, `nstressstrain`. The components of the stress vector are
  calculated and stored in the `stress` vector.
`output` =  array which is (if necessary) allocated  in an appropriate size, filled
  with the output quantity, and returned.

Note on the principal stresses: The principal stresses are calculated for the
fully three-dimensional  stress state, that is not the "in-plane" maximum and
minimum,  but rather  the  three-dimensional maximum (1) and minimum (3).
The intermediate principal stress is (2).
"""
function update2dstrn!(self::MatDeforElastOrtho,  stress::FFltVec, output::FFltVec, strain::FFltVec, thstrain::FFltVec=zeros(4), t::FFlt= 0.0, dt::FFlt= 0.0, loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)
    @assert length(stress) == nstressstrain(self.mr)
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
        copyto!(output, stress); output[4] = sz
    elseif quantity == :pressure || quantity == :Pressure
        (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
        sz = dot(self.D[3, 1:2], strain[1:2]-thstrain[1:2])-self.D[3,3]*thstrain[4];
        output[1] = -(sum(stress[[1,2]]) + sz)/3.
    elseif quantity == :princCauchy || quantity == :princcauchy
        (length(output) >= 3) || (output = zeros(3)) # make sure we can store it
        t = zeros(FFlt, 3,3)
        sz = dot(self.D[3, 1:2], strain[1:2]-thstrain[1:2])-self.D[3,3]*thstrain[4];
        t = stress4vto3x3t!(t, vcat(stress[1:3], [sz]));
        ep = eig(t);
        (length(output) >= 3) || (output = zeros(3)) # make sure we can store it
        copyto!(output,  sort(ep[1], rev=true));
    elseif quantity==:vonMises || quantity==:vonmises || quantity==:von_mises || quantity==:vm
        (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
        sz = dot(self.D[3, 1:2], strain[1:2]-thstrain[1:2])-self.D[3,3]*thstrain[4];
        s1=stress[1]; s2=stress[2]; s3=sz;
        s4=stress[3]; s5=0.0; s6=0.0;
        (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
        output[1] = sqrt(1.0/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)))
    end
    return output
end

"""
    thermalstrain2dstrn!(self::MatDeforElastOrtho, thstrain::FFltVec, dT= 0.0)

Compute thermal strain from the supplied temperature increment.

The thermal strain is evaluated  for the  two normal strains, the shear
strain, and finally for the through the thickness strain.
thstrain = [ex, ey, 0.0, ez].
"""
function thermalstrain2dstrn!(self::MatDeforElastOrtho, thstrain::FFltVec, dT= 0.0)
    @assert length(thstrain) == nthermstrain(self.mr)
    thstrain[1] = self.CTE1*dT
    thstrain[2] = self.CTE2*dT
    thstrain[3] = 0.0
    thstrain[4] = self.CTE3*dT
    return thstrain
end

################################################################################
# 2-D axially symmetric
################################################################################

"""
    tangentmoduli2daxi!(self::MatDeforElastOrtho,
      D::FFltMat,
      t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)

Calculate the material stiffness matrix.

`D` = matrix of tangent moduli, 3 x 3, supplied as a buffer and overwritten.
"""
function tangentmoduli2daxi!(self::MatDeforElastOrtho, D::FFltMat, t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)
    for i = 1:4
        for j = 1:4
            D[i, j] = self.D[i, j];
        end
    end
    return D
end

"""
    update2daxi!(self::MatDeforElastOrtho,  stress::FFltVec, output::FFltVec,
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
  stress and strain components, `nstressstrain`. The components of the stress vector are
  calculated and stored in the `stress` vector.
`output` =  array which is (if necessary) allocated  in an appropriate size, filled
  with the output quantity, and returned.
"""
function update2daxi!(self::MatDeforElastOrtho,  stress::FFltVec, output::FFltVec, strain::FFltVec, thstrain::FFltVec=zeros(3), t::FFlt= 0.0, dt::FFlt= 0.0, loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)
    @assert length(stress) == nstressstrain(self.mr)
    D = zeros(4, 4)
    tangentmoduli2daxi!(self, D, t, dt, loc, label)
    A_mul_B!(stress, D, strain-thstrain);
    if quantity == :nothing
        #Nothing to be copied to the output array
    elseif quantity == :cauchy || quantity == :Cauchy
        (length(output) >= 4) || (output = zeros(4)) # make sure we can store it
        copyto!(output, stress);
    elseif quantity == :pressure || quantity == :Pressure
        (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
        output[1] = -sum(stress[[1,2,3]])/3.
    elseif quantity == :princCauchy || quantity == :princcauchy
        t = zeros(FFlt,3,3)
        t = stress4vto3x3t!(t, stress[[1,2,4,3]]);
        ep = eig(t);
        (length(output) >= 3) || (output = zeros(3)) # make sure we can store it
        copyto!(output,  sort(ep[1], rev=true));
    elseif quantity==:vonMises || quantity==:vonmises || quantity==:von_mises || quantity==:vm
        s1=stress[1]; s2=stress[2]; s3=stress[3];
        s4=stress[4]; s5=0.0; s6=0.0;
        (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
        output[1] = sqrt(1.0/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)))
    end
    return output
end

"""
    thermalstrain2daxi!(self::MatDeforElastOrtho, thstrain::FFltVec, dT= 0.0)

Compute thermal strain from the supplied temperature increment.

The thermal strain is evaluated  for the  three normal strains and the shear
strain.
"""
function thermalstrain2daxi!(self::MatDeforElastOrtho, thstrain::FFltVec, dT= 0.0)
    @assert length(thstrain) == nthermstrain(self.mr)
    thstrain[1] = self.CTE1*dT
    thstrain[2] = self.CTE2*dT
    thstrain[3] = self.CTE3*dT
    thstrain[4] = 0.0
    return thstrain
end

end
