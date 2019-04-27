"""
    FEMMDeforLinearModule

Module for operations on interiors of domains to construct system matrices and
system vectors for linear deformation models.
"""
module FEMMDeforLinearModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FENodeSetModule: FENodeSet
import FinEtools.FESetModule: AbstractFESet, manifdim
import FinEtools.IntegDomainModule: IntegDomain
import FinEtools.FEMMDeforLinearBaseModule: AbstractFEMMDeforLinear
import FinEtools.DeforModelRedModule: AbstractDeforModelRed, DeforModelRed2DAxisymm
import FinEtools.MatDeforLinearElasticModule: AbstractMatDeforLinearElastic
import FinEtools.CSysModule: CSys

"""
    FEMMDeforLinear{MR<:AbstractDeforModelRed,  S<:AbstractFESet, F<:Function, M<:AbstractMatDeforLinearElastic} <: AbstractFEMMDeforLinear

Class for linear deformation finite element modeling machine.
"""
mutable struct FEMMDeforLinear{MR<:AbstractDeforModelRed,  S<:AbstractFESet, F<:Function, M<:AbstractMatDeforLinearElastic} <: AbstractFEMMDeforLinear
    mr::Type{MR} # model reduction type
    integdomain::IntegDomain{S, F} # integration domain data 
    mcsys::CSys # updater of the material orientation matrix
    material::M # material object
end
"""
    FEMMDeforLinear(mr::Type{MR}, integdomain::IntegDomain{S, F}, material::M) where {MR<:AbstractDeforModelRed, S<:AbstractFESet, F<:Function, M<:AbstractMatDeforLinearElastic}

Constructor of linear deformation finite element modeling machine.
"""
function FEMMDeforLinear(mr::Type{MR}, integdomain::IntegDomain{S, F}, material::M) where {MR<:AbstractDeforModelRed, S<:AbstractFESet, F<:Function, M<:AbstractMatDeforLinearElastic}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (integdomain.axisymmetric) || (mr != DeforModelRed2DAxisymm) "Axially symmetric requires axisymmetric to be true"
    return FEMMDeforLinear(mr, integdomain, CSys(manifdim(integdomain.fes)), material)
end

end
