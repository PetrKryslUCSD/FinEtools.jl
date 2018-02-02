"""
    FTypesModule

Module to define  basic types.

These are the types of quantities used in FinEtools. They are all exported.

export FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec
export FDataDict
"""
module FTypesModule

import Base.Complex

const FInt = Int
const FFlt = Float64
const FCplxFlt = Complex{FFlt}
FMat{T<:Number} = Array{T, 2}
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
