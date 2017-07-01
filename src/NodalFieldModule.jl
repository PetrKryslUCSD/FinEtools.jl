module NodalFieldModule

using FinEtools
using FinEtools.FTypesModule
using FinEtools.FieldModule
using FinEtools.FieldModule.@add_Field_fields

"""
    NodalField{T<:Number}

Nodal field.
"""
mutable struct NodalField{T<:Number} <: Field
  @add_Field_fields()
end
export NodalField

# Constructor of nodal field
function NodalField{T<:Number}(data::FMat{T}=[])
  values = deepcopy(data)
  dofnums = 0*similar(values,FInt)
  is_fixed = similar(values,Bool)
  fill!(is_fixed, 0)
  fixed_values = zeros(T,size(values))
  nfreedofs = 0
  return NodalField(values, dofnums, is_fixed, fixed_values, nfreedofs)
end

"""
    nnodes(self::NodalField)::FInt = nents(self)

Provide the number of nodes  in the nodal field.
"""
nnodes(self::NodalField)::FInt = nents(self)
export nnodes

end
