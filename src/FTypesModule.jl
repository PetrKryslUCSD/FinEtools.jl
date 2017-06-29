"""
Define  basic types.
"""
module FTypesModule

import Base.Complex

const FInt = Int
const FFlt = Float64
const FCplxFlt = Complex{Float64}
const FFltVec = Vector{FFlt}
const FIntVec = Vector{FInt}
const FFltMat = Matrix{FFlt}
const FIntMat = Matrix{FInt}
FMat{T<:Number} = Matrix{T}
FVec{T<:Number} = Vector{T}
export FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec

end
