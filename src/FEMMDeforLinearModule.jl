"""
    FEMMDeforLinearModule

Module for operations on interiors of domains to construct system matrices and
system vectors for linear deformation models.
"""
module FEMMDeforLinearModule

export FEMMDeforLinear
export stiffness, nzebcloadsstiffness, thermalstrainloads,
       inspectintegpoints

using FinEtools.FTypesModule
using FinEtools.FESetModule
using FinEtools.CSysModule
using FinEtools.GeoDModule
using FinEtools.FEMMBaseModule
using FinEtools.FEMMDeforLinearBaseModule
using FinEtools.FieldModule
using FinEtools.NodalFieldModule
using FinEtools.ElementalFieldModule
using FinEtools.ForceIntensityModule
using FinEtools.AssemblyModule
using FinEtools.DeforModelRedModule
using FinEtools.MatDeforModule
using FinEtools.FENodeToFEMapModule

"""
    FEMMDeforLinear{S<:FESet, F<:Function, P<:PropertyDeformationLinear}

Class for linear deformation finite element modeling machine.
"""
mutable struct FEMMDeforLinear{MR<:DeforModelRed,
  S<:FESet, F<:Function, M<:MatDefor} <: FEMMDeforLinearAbstract
  mr::Type{MR}
  geod::GeoD{S, F} # geometry data finite element modeling machine
  material::M # material object
  function FEMMDeforLinear(mr::Type{MR}, geod::GeoD{S, F},
    material::M) where {MR<:DeforModelRed,
    S<:FESet, F<:Function, M<:MatDefor}
    @assert mr === material.mr "Model reduction is mismatched"
    @assert (geod.axisymmetric) || (mr != DeforModelRed2DAxisymm) "Axially symmetric requires axisymmetric to be true"
    return new{MR, S, F, M}(mr, geod,  material)
  end
end


end
