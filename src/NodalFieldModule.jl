"""
    NodalFieldModule

Module for nodal fields.
"""
module NodalFieldModule

__precompile__(true)

import ..FieldModule: AbstractField, nents
import ..FieldModule.@add_Field_fields

"""
    NodalField{T<:Number, IT<:Integer} <: AbstractField

Nodal field, meaning the entities are the finite element nodes.
"""
mutable struct NodalField{T<:Number, IT<:Integer} <: AbstractField
    @add_Field_fields()
end

"""
    NodalField(data::Matrix{T}, zi::IT) where {T<:Number, IT<:Integer}

Constructor of nodal field. The values of the field are given by the array
on input, `data`. This array needs to have as many rows as there are nodes,
and as many columns as there are degrees of freedom per node.

The integer type for the storage of the degree of freedom numbers is set as that
of the argument `zi`.
"""
function NodalField(data::Matrix{T}, zi::IT) where {T<:Number, IT<:Integer}
    values = deepcopy(data)
    dofnums = 0 * similar(values, IT)
    is_fixed = similar(values, Bool)
    fill!(is_fixed, false)
    fixed_values = zeros(T, size(values))
    nfreedofs = zero(IT)
    return NodalField(values, dofnums, is_fixed, fixed_values, nfreedofs)
end

function NodalField(data::Matrix{T}) where {T<:Number}
    return NodalField(data, zero(Int))
end

"""
    NodalField(data::Vector{T}) where {T<:Number}

Constructor of nodal field. The values of the field are given by the vector on
input, `data`. This vector needs to have as many entries as there are nodes;
there is just one degree of freedom per node.
"""
function NodalField(data::Vector{T}) where {T<:Number}
    return NodalField(reshape(data, length(data), 1))
end

"""
    nnodes(self::NodalField)

Provide the number of nodes  in the nodal field.
"""
nnodes(self::NodalField) = nents(self)

end
