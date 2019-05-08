"""
    ElementalFieldModule

Module for elemental fields.
"""
module ElementalFieldModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FieldModule: AbstractField, nents
import FinEtools.FieldModule.@add_Field_fields


"""
    ElementalField{T<:Number} <: AbstractField

Elemental field, meaning the entities are finite elements.

The values in the field are indexed by the element number.  This means  that
there needs to be one field per finite element set.
"""
mutable struct ElementalField{T<:Number} <: AbstractField
	@add_Field_fields()
end

"""
   ElementalField(data::FMat{T}=[]) where {T<:Number}

Constructor of elemental field. The values of the field are given by the array
on input, `data`. This array needs to have as many rows as there are elements,
and as many columns as there are degrees of freedom per element.
"""
function ElementalField(data::FMat{T}=[]) where {T<:Number}
	values = deepcopy(data)
	dofnums = 0*similar(values,FInt)
	is_fixed = similar(values,Bool)
	fill!(is_fixed, 0)
	fixed_values = zeros(T,size(values))
	nfreedofs = 0
	return ElementalField(values, dofnums, is_fixed, fixed_values, nfreedofs)
end

"""
    ElementalField(data::FVec{T}) where {T<:Number}

Constructor of elemental field. The values of the field are given by the vector
on input, `data`. This vector needs to have as many entries as there are elements;
there is just one degree of freedom per element.
"""
function ElementalField(data::FVec{T}) where {T<:Number}
    return ElementalField(reshape(data, length(data), 1))
end 

"""
    nelems(self::ElementalField)::FInt = nents(self)

Provide the number of elements  in the elemental field.
"""
nelems(self::ElementalField)::FInt = nents(self)

end
