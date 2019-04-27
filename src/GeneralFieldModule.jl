"""
    GeneralFieldModule

Module for general fields.
"""
module GeneralFieldModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FieldModule.AbstractField
import FinEtools.FieldModule.@add_Field_fields

"""
    GeneralField{T<:Number} <: AbstractField

General field, meaning the entities can be anything.
"""
mutable struct GeneralField{T<:Number} <: AbstractField
    @add_Field_fields()
end


"""
    GeneralField(data::FMat{T}=[]) where {T<:Number}

Constructor of general field.  The values of the field are given by the array
on input, `data`. This array needs to have as many rows as there are entities,
and as many columns as there are degrees of freedom per entities.
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

"""
    GeneralField(data::FVec{T}) where {T<:Number}

Constructor of general field.  The values of the field are given by the vector
on input, `data`. This vector needs to have as many rows as there are entities.
"""
function GeneralField(data::FVec{T}) where {T<:Number}
    return GeneralField(reshape(data, length(data), 1))
end 

end
