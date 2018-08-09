module mocylpull14nnn
using FinEtools
using Test
using InteractiveUtils
function test()
    E1=1.0;
    E2=1.0;
    nu23=0.19;
    p= 0.15;
    rin=1.;
    rex =1.2;
    Length = 1*rex

    MR = DeforModelRed2DAxisymm
    fens,fes = Q4block(rex-rin,Length,5,20);
    material = MatDeforElastIso(MR, 00.0, E1, nu23, 0.0)
    
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material, true)
    println("========== With printing ==========")
    @show @code_lowered FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material, true)
    println("========== Original ==========")
    @show @code_lowered FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)
    
    true
end
end
using .mocylpull14nnn
mocylpull14nnn.test()


# module mocylpull14a
# using FinEtools
# using Test
# function test()
#     MR = DeforModelRed2DAxisymm
#     material = MatDeforElastIso(MR, 0.0, 1.0, 0.0, 0.0)
#     @test MR === material.mr
#     fens,fes = Q4block(0.2,1.2,5,20);
#     println("Success? This succeeds both locally and in CI")
#     # @code_llvm FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material, true)
#     femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material, true)
#     println("Failure? This succeeds locally, but fails in CI")
#     # @code_llvm FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)
#     femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)
# end
# end
# using .mocylpull14a
# mocylpull14a.test()
