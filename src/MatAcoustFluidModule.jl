"""
    MatAcoustFluidModule

Module for acoustic-fluid  material.
"""
module MatAcoustFluidModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.MatModule: AbstractMat

"""
    MatAcoustFluid <: AbstractMat

Type for acoustic fluid material.
"""
struct MatAcoustFluid <: AbstractMat
	bulk_modulus::FFlt;# Bulk modulus
	mass_density::FFlt;# Mass density
end

"""
    bulkmodulus(self::MatAcoustFluid)

Return the bulk modulus.
"""
function bulkmodulus(self::MatAcoustFluid)
	return self.bulk_modulus
end

end
