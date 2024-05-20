"""
    GeneralFieldModule

Module for general fields.
"""
module GeneralFieldModule

__precompile__(true)

import ..FieldModule.AbstractField
import ..FieldModule.@add_Field_fields
import ..FieldModule.DOF_KIND_FREE

"""
    GeneralField{T<:Number, IT<:Integer} <: AbstractField

General field, meaning the entities can be anything.
"""
mutable struct GeneralField{T<:Number,IT<:Integer} <: AbstractField
    @add_Field_fields()
end

"""
    GeneralField(data::Matrix{T}, zi::IT) where {T<:Number, IT<:Integer}

Constructor of general field.  The values of the field are given by the array on
input, `data`. This array needs to have as many rows as there are entities, and
as many columns as there are degrees of freedom per entities.

The integer type for the storage of the degree of freedom numbers is set as that
of the argument `zi`.
"""
function GeneralField(data::Matrix{T}, zi::IT) where {T<:Number,IT<:Integer}
    values = deepcopy(data)
    dofnums = 0 * similar(values, IT)
    kind = similar(values, Int8)
    fill!(kind, DOF_KIND_FREE)
    return GeneralField(values, dofnums, kind, UnitRange{IT}[])
end

function GeneralField(data::Matrix{T}) where {T<:Number}
    return GeneralField(data, zero(Int))
end

"""
    GeneralField(data::Vector{T}) where {T<:Number}

Constructor of general field.  The values of the field are given by the vector
on input, `data`. This vector needs to have as many rows as there are
entities.
"""
function GeneralField(data::Vector{T}) where {T<:Number}
    return GeneralField(reshape(data, length(data), 1))
end

end
