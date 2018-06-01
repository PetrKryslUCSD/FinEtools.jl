"""
    GeneralFieldModule

Module for general fields.
"""
module GeneralFieldModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FieldModule.Field
import FinEtools.FieldModule.@add_Field_fields

"""
    GeneralField{T<:Number} <: Field

General field.
"""
mutable struct GeneralField{T<:Number} <: Field
    @add_Field_fields()
end


"""
    GeneralField(data::FMat{T}=[]) where {T<:Number}

Constructor of general field.
"""
function GeneralField(data::FMat{T}=[]) where {T<:Number}
    values = deepcopy(data)
    dofnums = 0*similar(values,FInt)
    is_fixed = similar(values,Bool)
    fill!(is_fixed, 0)
    fixed_values = zeros(T,size(values))
    nfreedofs = 0
    return GeneralField(values, dofnums, is_fixed, fixed_values, nfreedofs)
end

function GeneralField(data::FVec{T}) where {T<:Number}
    return GeneralField(reshape(data, length(data), 1))
end 

end
