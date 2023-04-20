"""
    MatModule

Module for abstract material.
"""
module MatModule

__precompile__(true)

"""
    AbstractMat

Abstract type of material.
"""
abstract type AbstractMat end

"""
    massdensity(self::AbstractMat)

Return mass density.
"""
function massdensity(self::AbstractMat)
    return self.mass_density
end

end
