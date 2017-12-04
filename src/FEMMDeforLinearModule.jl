"""
    FEMMDeforLinearModule

Module for operations on interiors of domains to construct system matrices and
system vectors for linear deformation models.
"""
module FEMMDeforLinearModule

# export FEMMDeforLinear
# export stiffness, nzebcloadsstiffness, thermalstrainloads,
#        inspectintegpoints

using FinEtools

"""
    FEMMDeforLinear{S<:FESet, F<:Function, P<:PropertyDeformationLinear}

Class for linear deformation finite element modeling machine.
"""
mutable struct FEMMDeforLinear{MR<:DeforModelRed,  S<:FESet, F<:Function, M<:MatDefor} <: FEMMDeforLinearAbstract
    mr::Type{MR}
    integdata::IntegData{S, F} # geometry data 
    mcsys::CSys # updater of the material orientation matrix
    material::M # material object
end

function FEMMDeforLinear(mr::Type{MR}, integdata::IntegData{S, F}, material::M) where {MR<:DeforModelRed, S<:FESet, F<:Function, M<:MatDefor}
    @assert mr === material.mr "Model reduction is mismatched"
    @assert (integdata.axisymmetric) || (mr != DeforModelRed2DAxisymm) "Axially symmetric requires axisymmetric to be true"
    return FEMMDeforLinear(mr, integdata, CSys(manifdim(integdata.fes)), material)
end

end
