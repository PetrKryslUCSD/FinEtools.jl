"""
    FEMMDeforLinearModule

Module for operations on interiors of domains to construct system matrices and
system vectors for linear deformation models.
"""
module FEMMDeforLinearModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FENodeSetModule: FENodeSet
import FinEtools.FESetModule: FESet, manifdim
import FinEtools.IntegDomainModule: IntegDomain
import FinEtools.FEMMDeforLinearBaseModule: FEMMDeforLinearAbstract
import FinEtools.DeforModelRedModule: DeforModelRed, DeforModelRed2DAxisymm
import FinEtools.MatDeforModule: MatDefor
import FinEtools.CSysModule: CSys

"""
    FEMMDeforLinear{S<:FESet, F<:Function, P<:PropertyDeformationLinear}

Class for linear deformation finite element modeling machine.
"""
mutable struct FEMMDeforLinear{MR<:DeforModelRed,  S<:FESet, F<:Function, M<:MatDefor} <: FEMMDeforLinearAbstract
    mr::Type{MR}
    integdomain::IntegDomain{S, F} # geometry data 
    mcsys::CSys # updater of the material orientation matrix
    material::M # material object
end

function FEMMDeforLinear(mr::Type{MR}, integdomain::IntegDomain{S, F}, material::M) where {MR<:DeforModelRed, S<:FESet, F<:Function, M<:MatDefor}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (integdomain.axisymmetric) || (mr != DeforModelRed2DAxisymm) "Axially symmetric requires axisymmetric to be true"
    return FEMMDeforLinear(mr, integdomain, CSys(manifdim(integdomain.fes)), material)
end

end
