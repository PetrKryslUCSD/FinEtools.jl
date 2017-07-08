"""
    MatHeatDiffModule

Module for linear heat diffusion material models.
"""
module MatHeatDiffModule

export MatHeatDiff

using FinEtools.FTypesModule

# Type for heat diffusion models of materials.
struct MatHeatDiff
  thermal_conductivity::FFltMat;# Thermal conductivity
  specific_heat::FFlt;# Specific heat per unit volume
end


function MatHeatDiff(thermal_conductivity)
    return MatHeatDiff(thermal_conductivity, 0.0)
end

# """
#     bar
#
# Compute.
# """
# function update!(self, ms, gradtheta, Rm; output = nothing, outputRm =nothing)
#
#   kappa = self.property.thermal_conductivity;
#   gradtheta = reshape(context.gradtheta, [],1);# make column vector
#   flux = - kappa * gradtheta;# in local material orientation coordinates
#   if outputRm!= nothing #  output coordinate system supplied?
#     if Rm!=nothing
#       flux = (Rm*flux);# in global coordinate system
#     end
#     flux = outputRm'*flux;# in output coordinate system
#   end
#   if output!= nothing
#     if output=="flux"
#       out = flux;
#     else
#       out = nothing;
#     end
#   else
#     out = flux;
#   end
#   return out, ms
# end


end
