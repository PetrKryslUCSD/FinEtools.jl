"""
    MatAcoustFluidModule

Module for acoustic-fluid  material.
"""
module MatAcoustFluidModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict

# Class for acoustic fluid models of Mats.
struct MatAcoustFluid
  bulk_modulus::FFlt;# Bulk modulus
  mass_density::FFlt;# Mass density
end


# function update(self, ms, gradtheta, Rm; output = nothing, outputRm =nothing)
#   # # _Update_ material state.
#   # #
#   # # function [out, newms] = update (self, ms, context)
#   # #
#   # # Update material state.  Return the updated material state, and the
#   # # requested quantity (default is the heat flux).
#   # #
#   # # Arguments
#   # #     m=material
#   # #     ms = material state
#   # #     gradtheta= temperature gradient in the local material
#   # #           directions (which may be the same as the global coordinate
#   # #           directions)
#   # #     Rm = material orientation  matrix, may be supplied as empty  if it
#   # #           corresponds to an identity.
#   # #     and optional arguments
#   # #         output=type of quantity to output, and interpreted by the
#   # #           particular material; [] is returned when the material does not
#   # #           recognize the requested quantity to indicate uninitialized
#   # #           value.  It can be tested with isempty ().
#   # #              output ='flux' - heat flux vector in the local material
#   # #                      directions; this is the default
#   # #                      when output type is not specified.
#   # #         outputRm = orientation matrix  of the coordinate system in which
#   # #           the output should be calculated;  if this matrix is not supplied,
#   # #           it is assumed that the output is to be provided  in the local
#   # #           material coordinate system; otherwise the output is first
#   # #           transformed to the global coordinate system, and then to the
#   # #           output coordinate system.
#   # # Output
#   # #     out=requested quantity
#   # #     newms=new material state; don't forget that if the update is final
#   # #           the material state newms must be assigned and stored.  Otherwise
#   # #           the material update is lost!
#   # #
#   kappa = self.property.thermal_conductivity;
#   gradtheta = reshape (context.gradtheta, [],1);# make column vector
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
