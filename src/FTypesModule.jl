"""
    FTypesModule

Module to define  basic types.
"""
module FTypesModule

import Base.Complex

const FInt = Int
const FFlt = Float64
const FCplxFlt = Complex{Float64}
FMat{T<:Number} = Matrix{T}
FVec{T<:Number} = Vector{T}
const FFltVec = FVec{FFlt}
const FIntVec = FVec{FInt}
const FFltMat = FMat{FFlt}
const FIntMat = FMat{FInt}
export FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec

end
