module ElementalFieldModule

using FinEtools
using FinEtools.FTypesModule
using FinEtools.FieldModule
using FinEtools.FieldModule.@add_Field_fields


"""
    ElementalField{T<:Number}

Elemental field.

The values in the field are indexed by the element number.  This means  that
there needs to be one field per finite element set.
"""
type ElementalField{T<:Number} <: Field
  @add_Field_fields()
end
export ElementalField

# Constructor of nodal field
function ElementalField{T<:Number}(data::FMat{T}=[])
  values = deepcopy(data)
  dofnums = 0*similar(values,FInt)
  is_fixed = similar(values,Bool)
  fill!(is_fixed, 0)
  fixed_values = zeros(T,size(values))
  nfreedofs = 0
  return ElementalField(values, dofnums, is_fixed, fixed_values, nfreedofs)
end

"""
    nelems(self::ElementalField)::FInt = nents(self)

Provide the number of elements  in the elemental field.
"""
nelems(self::ElementalField)::FInt = nents(self)
export nelems

end
