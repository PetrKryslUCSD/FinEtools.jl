"""
    ElementalFieldModule

Module for elemental fields.
"""
module ElementalFieldModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FieldModule: Field, nents
import FinEtools.FieldModule.@add_Field_fields


"""
    ElementalField{T<:Number}

Elemental field.

The values in the field are indexed by the element number.  This means  that
there needs to be one field per finite element set.
"""
mutable struct ElementalField{T<:Number} <: Field
  @add_Field_fields()
end

"""
   ElementalField(data::FMat{T}=[]) where {T<:Number}

Constructor of elemental field.
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

function ElementalField(data::FVec{T}) where {T<:Number}
    return ElementalField(reshape(data, length(data), 1))
end 

"""
    nelems(self::ElementalField)::FInt = nents(self)

Provide the number of elements  in the elemental field.
"""
nelems(self::ElementalField)::FInt = nents(self)

end
