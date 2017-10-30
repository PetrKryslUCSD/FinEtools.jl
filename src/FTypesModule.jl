"""
    FTypesModule

Module to define  basic types.
"""
module FTypesModule

export FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec
export FDataDict

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

"""
    FDataDict = Dict{String, Any}

Type for the model-data packaging system (used by all FinEtools algorithms).  
"""
const FDataDict = Dict{String, Any}

end
