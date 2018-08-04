module mocylpull14a
using FinEtools
using Test
function test()
    MR = DeforModelRed2DAxisymm
    material = MatDeforElastIso(MR, 0.0, 1.0, 0.0, 0.0)
    @test MR === material.mr
    fens,fes = Q4block(0.2,1.2,5,20);
    println("Success? This succeeds both locally and in CI")
    # @code_llvm FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material, true)
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material, true)
    println("Failure? This succeeds locally, but fails in CI")
    # @code_llvm FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)
end
end
using .mocylpull14a
mocylpull14a.test()