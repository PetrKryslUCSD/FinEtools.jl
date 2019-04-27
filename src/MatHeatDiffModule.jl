"""
    MatHeatDiffModule

Module for linear heat diffusion material models.
"""
module MatHeatDiffModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.MatModule: AbstractMat
import LinearAlgebra: mul!, Transpose
At_mul_B!(C, A, B) = mul!(C, Transpose(A), B)
A_mul_B!(C, A, B) = mul!(C, A, B)

"""
    MatHeatDiff{MTAN<:Function, MUPD<:Function} <: AbstractMat

Type of material model for heat diffusion.
"""
struct MatHeatDiff{MTAN<:Function, MUPD<:Function} <: AbstractMat
	thermal_conductivity::Array{FFlt, 2};# Thermal conductivity
	specific_heat::FFlt;# Specific heat per unit volume
	mass_density::FFlt # mass density
	tangentmoduli!::MTAN
	update!::MUPD
end

"""
    MatHeatDiff(thermal_conductivity)

Construct material model for heat diffusion.

Supply the matrix of thermal conductivity constants.
"""
function MatHeatDiff(thermal_conductivity)
    return MatHeatDiff(thermal_conductivity, 0.0, 0.0, tangentmoduli!, update!)
end

"""
    tangentmoduli!(self::MatHeatDiff, kappabar::FFltMat, t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)

Calculate the thermal conductivity matrix.

- `kappabar` = matrix of thermal conductivity (tangent moduli) in material
  coordinate system, supplied as a buffer and overwritten.
"""
function tangentmoduli!(self::MatHeatDiff, kappabar::FFltMat, t::FFlt = 0.0, dt::FFlt = 0.0, loc::FFltMat = reshape(FFlt[],0,0), label::FInt = 0)
    copyto!(kappabar, self.thermal_conductivity);
    return kappabar
end

"""
    update!(self::MatHeatDiff, heatflux::FFltVec, output::FFltVec, gradT::FFltVec, t::FFlt= 0.0, dt::FFlt= 0.0,  loc::FFltMat=FFltMat[], label::FInt=0, quantity=:nothing)

Update material state.

- `strain` = strain vector,
`thstrain` = thermal strain vector,
- `t` = current time,
- `dt` = current time step,
- `loc` = location of the quadrature point in global Cartesian coordinates,
- `label` = label of the finite element in which the quadrature point is located.

Output:
- `heatflux` = heat flux vector, allocated by the caller with a size of the
  embedding space. The components of the heat flux vector are calculated and
  stored in the `heatflux` vector.
- `output` =  array which is (if necessary) allocated  in an appropriate size,
  filled with the output quantity, and returned.
"""
function update!(self::MatHeatDiff, heatflux::FFltVec, output::FFltVec, gradT::FFltVec, t::FFlt= 0.0, dt::FFlt= 0.0, loc::FFltMat=reshape(FFlt[],0,0), label::FInt=0, quantity=:nothing)
	sdim = size(self.thermal_conductivity, 2)
    @assert length(heatflux) == sdim
    A_mul_B!(heatflux, self.thermal_conductivity, -gradT);
    if quantity == :nothing
        #Nothing to be copied to the output array
    elseif quantity == :heatflux 
        (length(output) >= sdim) || (output = zeros(sdim)) # make sure we can store it
        copyto!(output, heatflux);
    end
    return output
end

end
