"""
    ElementalFieldModule

Module for elemental fields.
"""
module ElementalFieldModule

__precompile__(true)

import ..FieldModule: AbstractField, nents
import ..FieldModule.@add_Field_fields
import ..FieldModule.KIND_INT
import ..FieldModule.DOF_KIND_FREE

"""
    ElementalField{T<:Number, IT<:Integer} <: AbstractField

Elemental field, meaning the entities are finite elements.

The values in the field are indexed by the element number.  This means  that
there needs to be one field per finite element set.
"""
mutable struct ElementalField{T<:Number,IT<:Integer} <: AbstractField
    @add_Field_fields()
end

"""
   ElementalField(data::Matrix{T}, zi::IT) where {T<:Number, IT<:Integer}

Constructor of elemental field. The values of the field are given by the array
on input, `data`. This array needs to have as many rows as there are elements,
and as many columns as there are degrees of freedom per element.

The integer type for the storage of the degree of freedom numbers is set as that
of the argument `zi`.
"""
function ElementalField(data::Matrix{T}, zi::IT) where {T<:Number,IT<:Integer}
    values = deepcopy(data)
    dofnums = 0 * similar(values, IT)
    kind = similar(values, Int8)
    fill!(kind, DOF_KIND_FREE)
    return ElementalField(values, dofnums, kind, UnitRange{IT}[])
end

function ElementalField(data::Matrix{T}) where {T<:Number}
    return ElementalField(data, zero(Int))
end

"""
    ElementalField(data::Vector{T}) where {T<:Number}

Constructor of elemental field. The values of the field are given by the vector
on input, `data`. This vector needs to have as many entries as there are elements;
there is just one degree of freedom per element.
"""
function ElementalField(data::Vector{T}) where {T<:Number}
    return ElementalField(reshape(data, length(data), 1))
end

"""
    nelems(self::ElementalField)

Provide the number of elements  in the elemental field.
"""
nelems(self::ElementalField) = nents(self)

end
