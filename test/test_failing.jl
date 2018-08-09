module mocylpull14nnn
using FinEtools
using Test
using InteractiveUtils
function test()
    E1=1.0;
    nu23=0.19;
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