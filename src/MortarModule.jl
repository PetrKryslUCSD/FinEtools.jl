module MortarModule

__precompile__(true)

using LinearAlgebra
using SparseArrays

using ..FENodeSetModule: FENodeSet
using ..FESetModule: count
using ..IntegRuleModule: TriRule, GaussRule
using ..IntegDomainModule: IntegDomain
using ..NodalFieldModule: NodalField
using ..ElementalFieldModule: ElementalField
using ..FEMMBaseModule: FEMMBase, bilform_masslike

export common_refinement,
       clip_polygon,
       build_common_refinement_coupling

function common_refinement()
    # placeholder
    println("common_refinement is not implemented yet.")
    return nothing
end

end