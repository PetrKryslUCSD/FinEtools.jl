"""
    NodalFieldModule

Module for nodal fields.
"""
module NodalFieldModule

using ..FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import ..FieldModule: AbstractField, nents
import ..FieldModule.@add_Field_fields

"""
    NodalField{T<:Number} <: AbstractField

Nodal field, meaning the entities are the finite element nodes.
"""
mutable struct NodalField{T<:Number} <: AbstractField
    @add_Field_fields()
end

"""
    NodalField(data::FMat{T}=[]) where {T<:Number}

Constructor of nodal field. The values of the field are given by the array
on input, `data`. This array needs to have as many rows as there are nodes,
and as many columns as there are degrees of freedom per node.
"""
function NodalField(data::FMat{T}=[]) where {T<:Number}
    values = deepcopy(data)
    dofnums = 0*similar(values,FInt)
    is_fixed = similar(values,Bool)
    fill!(is_fixed, 0)
    fixed_values = zeros(T,size(values))
    nfreedofs = 0
    return NodalField(values, dofnums, is_fixed, fixed_values, nfreedofs)
end

"""
    NodalField(data::FVec{T}) where {T<:Number}

Constructor of nodal field. The values of the field are given by the vector
on input, `data`. This vector needs to have as many entries as there are nodes;
there is just one degree of freedom per nodes.
"""
function NodalField(data::FVec{T}) where {T<:Number}
    return NodalField(reshape(data, length(data), 1))
end 

"""
    nnodes(self::NodalField)::FInt = nents(self)

Provide the number of nodes  in the nodal field.
"""
nnodes(self::NodalField)::FInt = nents(self)

end
